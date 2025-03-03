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
! !ROUTINE: FLake1_main
! \label{FLake1_main}
!
! !REVISION HISTORY:
!
!   6/4/13: Shugong Wang; initial implementation for FLake1 with LIS-7
!
! !INTERFACE:
subroutine FLake1_main(n)
! !USES:
  use LIS_coreMod
  use LIS_histDataMod
  use LIS_timeMgrMod, only : LIS_isAlarmRinging
  use LIS_logMod, only     : LIS_logunit, LIS_endrun
  use FLake1_Mod
   !use other modules
  
  implicit none
! !ARGUMENTS:
  integer, intent(in)  :: n
  integer              :: t
  integer              :: i
  real                 :: dt
  real                 :: lat, lon
  integer              :: row, col
  logical              :: alarmCheck

!
! !DESCRIPTION:
!  This is the entry point for calling the FLake1 physics.
!  This routine calls the {\tt LIS7_FLake_2013 } routine that performs the
!  land surface computations, to solve for water and energy equations.

!  The arguments are:
!  \begin{description}
!  \item[n]
!    index of the nest
!  \end{description}
!EOP

! define variables for FLake1
  real             :: tmp_swdown             ! incoming solar radiation at the surface [W/m2]
  real             :: tmp_lwdown             ! incoming longwave radiation at the surface [W/m2]
  real             :: tmp_wind_e             ! eastward wind speed at height_wind [m/s]
  real             :: tmp_wind_n             ! northward wind speed at height_wind [m/s]
  real             :: tmp_tair               ! air temperature at height_tq [K]
  real             :: tmp_qair               ! specific air humidity at height_tq [kg/kg]
  real             :: tmp_psurf              ! surface air pressure [Pa]
  real             :: tmp_height_wind        ! height where wind speed is measured [m]
  real             :: tmp_height_tq          ! height where temperature and humidity are measured [-]
  real             :: tmp_flake_dt           ! model time step [s]
  real             :: tmp_lon                ! longitude of lake center [-]
  real             :: tmp_lat                ! latitude of lake center [-]
  real             :: tmp_depth_w            ! lake depth [m]
  real             :: tmp_fetch              ! typical wind fetch [m]
  real             :: tmp_depth_bs           ! depth of the thermally active layer of the bottom sediments [m]
  real             :: tmp_T_bs               ! temperature at the outer edge of the thermally active layer of the bottom sediments [K]
  real             :: tmp_T_snow             ! temperature at the air-snow interface [K]
  real             :: tmp_T_ice              ! temperature at the snow-ice interface [K]
  real             :: tmp_T_mnw              ! mean temperature of the water column [K]
  real             :: tmp_T_wML              ! temperature of mixed layer [K]
  real             :: tmp_T_bot              ! temperature at the water-bottom sediment interface [K]
  real             :: tmp_T_b1               ! temperature at the bottom of the upper layer of the sediments [K]
  real             :: tmp_C_T                ! thermocline shape factor [-]
  real             :: tmp_H_snow             ! snow thickness [m]
  real             :: tmp_H_ice              ! ice thickness [m]
  real             :: tmp_H_ML               ! thickness of mixed layer [m]
  real             :: tmp_H_B1               ! thickness of the upper layer of bottom sediments [m]
  real             :: tmp_T_sfc              ! surface temperature [K]
  real             :: tmp_albedo_water       ! water surface albedo with resect to solar radiation [-]
  real             :: tmp_albedo_ice         ! ice surface albedo with respect to the solar radiation [-]
  real             :: tmp_albedo_snow        ! snow surface albedo with respect to the solar radiation [-]
  real             :: tmp_ufr_a              ! friction velocity in air [m/s]
  real             :: tmp_ufr_w              ! friction velocity in surface water [m/s]
  real             :: tmp_Wconv              ! convective velocity scale [m/s]
  real             :: tmp_Q_se               ! sensible surface heat flux [W/m2]
  real             :: tmp_Q_la               ! latent surface heat flux [W/m2]
  real             :: tmp_I_w                ! radiation flux through the ice-water or air-water interface [W/m2]
  real             :: tmp_Q_lwa              ! longwave radiation flux from atmosphere [W/m2]
  real             :: tmp_Q_lww              ! longwave radiation flux from water [W/m2]
  real             :: tmp_Q_bot              ! heat flux across water-sediments boundary [W/m2]
  real             :: topog, temp, temp1
  character*3   :: fnest

  write(fnest,'(i3.3)') n

    ! check FLake1 alarm. If alarm is ring, run model. 
  alarmCheck = LIS_isAlarmRinging(LIS_rc, "FLAKE1 model alarm"//trim(fnest))
  if (alarmCheck) Then
     do t = 1, LIS_rc%npatch(n, LIS_rc%lake_index)
        

        dt = LIS_rc%ts
        row = LIS_surface(n, LIS_rc%lake_index)%tile(t)%row
        col = LIS_surface(n, LIS_rc%lake_index)%tile(t)%col
        lat = LIS_domain(n)%grid(LIS_domain(n)%gindex(col, row))%lat
        lon = LIS_domain(n)%grid(LIS_domain(n)%gindex(col, row))%lon

        ! retrieve forcing data from FLAKE1_struc(n)%flake1(t) and assign to local variables
        tmp_swdown = FLAKE1_struc(n)%flake1(t)%swdown / FLAKE1_struc(n)%forc_count
        tmp_lwdown = FLAKE1_struc(n)%flake1(t)%lwdown / FLAKE1_struc(n)%forc_count
        tmp_wind_e = FLAKE1_struc(n)%flake1(t)%wind_e / FLAKE1_struc(n)%forc_count
        tmp_wind_n = FLAKE1_struc(n)%flake1(t)%wind_n / FLAKE1_struc(n)%forc_count
        tmp_tair   = FLAKE1_struc(n)%flake1(t)%tair   / FLAKE1_struc(n)%forc_count
        tmp_qair   = FLAKE1_struc(n)%flake1(t)%qair   / FLAKE1_struc(n)%forc_count
        tmp_psurf  = FLAKE1_struc(n)%flake1(t)%psurf  / FLAKE1_struc(n)%forc_count

        if (trim(LIS_rc%startcode) .eq. "coldstart".and.&
             FLAKE1_struc(n)%startFlag ) then
           FLAKE1_struc(n)%startFlag = .false. 
           FLAKE1_struc(n)%flake1(t)%T_wML = tmp_tair
           FLAKE1_struc(n)%flake1(t)%T_mnw = tmp_tair
           FLAKE1_struc(n)%flake1(t)%T_bot = tmp_tair
           if(lat.lt.30) then 
              tmp_T_bot =  tmp_tair
           else
              topog = 0.0 ! for now
              tmp_t_bot = 28.9 - 0.43*lat - 0.0038*topog + 273.15
           endif

        endif
        ! check validity of swdown
        if(tmp_swdown .eq. LIS_rc%udef) then
           write(LIS_logunit, *) "undefined value found for forcing variable swdown in FLake1"
           write(LIS_logunit, *) "for tile ", t, "latitude = ", lat, "longitude = ", lon
           call LIS_endrun()
        endif
        ! check validity of lwdown
        if(tmp_lwdown .eq. LIS_rc%udef) then
           write(LIS_logunit, *) "undefined value found for forcing variable lwdown in FLake1"
           write(LIS_logunit, *) "for tile ", t, "latitude = ", lat, "longitude = ", lon
           call LIS_endrun()
        endif
        ! check validity of wind_e
        if(tmp_wind_e .eq. LIS_rc%udef) then
           write(LIS_logunit, *) "undefined value found for forcing variable wind_e in FLake1"
           write(LIS_logunit, *) "for tile ", t, "latitude = ", lat, "longitude = ", lon
           call LIS_endrun()
        endif
        ! check validity of wind_n
        if(tmp_wind_n .eq. LIS_rc%udef) then
           write(LIS_logunit, *) "undefined value found for forcing variable wind_n in FLake1"
           write(LIS_logunit, *) "for tile ", t, "latitude = ", lat, "longitude = ", lon
           call LIS_endrun()
        endif
        ! check validity of tair
        if(tmp_tair .eq. LIS_rc%udef) then
           write(LIS_logunit, *) "undefined value found for forcing variable tair in FLake1"
           write(LIS_logunit, *) "for tile ", t, "latitude = ", lat, "longitude = ", lon
           call LIS_endrun()
        endif
        ! check validity of qair
        if(tmp_qair .eq. LIS_rc%udef) then
           write(LIS_logunit, *) "undefined value found for forcing variable qair in FLake1"
           write(LIS_logunit, *) "for tile ", t, "latitude = ", lat, "longitude = ", lon
           call LIS_endrun()
        endif
        ! check validity of psurf
        if(tmp_psurf .eq. LIS_rc%udef) then
           write(LIS_logunit, *) "undefined value found for forcing variable psurf in FLake1"
           write(LIS_logunit, *) "for tile ", t, "latitude = ", lat, "longitude = ", lon
           call LIS_endrun()
        endif
        ! get parameters 
        tmp_height_wind    = FLAKE1_struc(n)%height_wind
        tmp_height_tq      = FLAKE1_struc(n)%height_tq  
        tmp_flake_dt       = FLAKE1_struc(n)%flake_dt   
        tmp_lon            = FLAKE1_struc(n)%flake1(t)%lon        
        tmp_lat            = FLAKE1_struc(n)%flake1(t)%lat        
        tmp_depth_w        = FLAKE1_struc(n)%flake1(t)%depth_w    
        if(tmp_depth_w.gt.60.0)  then !limit of flake
           tmp_depth_w = 60.0
        endif
        tmp_fetch          = FLAKE1_struc(n)%flake1(t)%fetch      
        tmp_depth_bs       = FLAKE1_struc(n)%flake1(t)%depth_bs   
        tmp_T_bs           = FLAKE1_struc(n)%flake1(t)%T_bs       
        ! get state variables
        tmp_T_snow          = FLAKE1_struc(n)%flake1(t)%T_snow      
        tmp_T_ice           = FLAKE1_struc(n)%flake1(t)%T_ice       
        tmp_T_mnw           = FLAKE1_struc(n)%flake1(t)%T_mnw       
        tmp_T_wML           = FLAKE1_struc(n)%flake1(t)%T_wML

        tmp_T_bot           = FLAKE1_struc(n)%flake1(t)%T_bot       

        tmp_T_b1            = FLAKE1_struc(n)%flake1(t)%T_b1        
        tmp_C_T             = FLAKE1_struc(n)%flake1(t)%C_T         
        tmp_H_snow          = FLAKE1_struc(n)%flake1(t)%H_snow      
        tmp_H_ice           = FLAKE1_struc(n)%flake1(t)%H_ice       
        tmp_H_ML            = FLAKE1_struc(n)%flake1(t)%H_ML        
        tmp_H_B1            = FLAKE1_struc(n)%flake1(t)%H_B1        
        tmp_T_sfc           = FLAKE1_struc(n)%flake1(t)%T_sfc       
        tmp_albedo_water    = FLAKE1_struc(n)%flake1(t)%albedo_water

        temp = 2*3.14159265*(LIS_rc%doy-1)/365.0
        temp1 = 0.006918-0.399912*cos(temp)+0.070257*sin(temp)-  &
             0.006758*cos(2.0*temp)+0.000907*sin(2.0*temp) -  &
             0.002697*cos(3.0*temp)+0.00148*sin(3.0*temp) 
        tmp_albedo_water = 0.06/cos((LAT*3.14159265/180.0-temp1)/1.2)

        tmp_albedo_ice      = FLAKE1_struc(n)%flake1(t)%albedo_ice  
        tmp_albedo_snow     = FLAKE1_struc(n)%flake1(t)%albedo_snow 

        ! call model physics 
        call LIS7_FLake_2013(t, tmp_swdown            , & ! in    - incoming solar radiation at the surface [W/m2]
             tmp_lwdown            , & ! in    - incoming longwave radiation at the surface [W/m2]
             tmp_wind_e            , & ! in    - eastward wind speed at height_wind [m/s]
             tmp_wind_n            , & ! in    - northward wind speed at height_wind [m/s]
             tmp_tair              , & ! in    - air temperature at height_tq [K]
             tmp_qair              , & ! in    - specific air humidity at height_tq [kg/kg]
             tmp_psurf             , & ! in    - surface air pressure [Pa]
             tmp_height_wind       , & ! in    - height where wind speed is measured [m]
             tmp_height_tq         , & ! in    - height where temperature and humidity are measured [-]
             tmp_flake_dt          , & ! in    - model time step [s]
             tmp_lon               , & ! in    - longitude of lake center [-]
             tmp_lat               , & ! in    - latitude of lake center [-]
             tmp_depth_w           , & ! in    - lake depth [m]
             tmp_fetch             , & ! in    - typical wind fetch [m]
             tmp_depth_bs          , & ! in    - depth of the thermally active layer of the bottom sediments [m]
             tmp_T_bs              , & ! in    - temperature at the outer edge of the thermally active layer of the bottom sediments [K]
             tmp_T_snow            , & ! inout - temperature at the air-snow interface [K]
             tmp_T_ice             , & ! inout - temperature at the snow-ice interface [K]
             tmp_T_mnw             , & ! inout - mean temperature of the water column [K]
             tmp_T_wML             , & ! inout - temperature of mixed layer [K]
             tmp_T_bot             , & ! inout - temperature at the water-bottom sediment interface [K]
             tmp_T_b1              , & ! inout - temperature at the bottom of the upper layer of the sediments [K]
             tmp_C_T               , & ! inout - thermocline shape factor [-]
             tmp_H_snow            , & ! inout - snow thickness [m]
             tmp_H_ice             , & ! inout - ice thickness [m]
             tmp_H_ML              , & ! inout - thickness of mixed layer [m]
             tmp_H_B1              , & ! inout - thickness of the upper layer of bottom sediments [m]
             tmp_T_sfc             , & ! inout - surface temperature [K]
             tmp_albedo_water      , & ! inout - water surface albedo with resect to solar radiation [-]
             tmp_albedo_ice        , & ! inout - ice surface albedo with respect to the solar radiation [-]
             tmp_albedo_snow       , & ! inout - snow surface albedo with respect to the solar radiation [-]
             tmp_ufr_a             , & ! out   - friction velocity in air [m/s]
             tmp_ufr_w             , & ! out   - friction velocity in surface water [m/s]
             tmp_Wconv             , & ! out   - convective velocity scale [m/s]
             tmp_Q_se              , & ! out   - sensible surface heat flux [W/m2]
             tmp_Q_la              , & ! out   - latent surface heat flux [W/m2]
             tmp_I_w               , & ! out   - radiation flux through the ice-water or air-water interface [W/m2]
             tmp_Q_lwa             , & ! out   - longwave radiation flux from atmosphere [W/m2]
             tmp_Q_lww             , & ! out   - longwave radiation flux from water [W/m2]
             tmp_Q_bot             )   ! out   - heat flux across water-sediments boundary [W/m2]

        ! save state variables from local variables to global variables
        FLAKE1_struc(n)%flake1(t)%T_snow          = tmp_T_snow         
        FLAKE1_struc(n)%flake1(t)%T_ice           = tmp_T_ice          
        FLAKE1_struc(n)%flake1(t)%T_mnw           = tmp_T_mnw          
        FLAKE1_struc(n)%flake1(t)%T_wML           = tmp_T_wML          
        FLAKE1_struc(n)%flake1(t)%T_bot           = tmp_T_bot          
        FLAKE1_struc(n)%flake1(t)%T_b1            = tmp_T_b1           
        FLAKE1_struc(n)%flake1(t)%C_T             = tmp_C_T            
        FLAKE1_struc(n)%flake1(t)%H_snow          = tmp_H_snow         
        FLAKE1_struc(n)%flake1(t)%H_ice           = tmp_H_ice          
        FLAKE1_struc(n)%flake1(t)%H_ML            = tmp_H_ML           
        FLAKE1_struc(n)%flake1(t)%H_B1            = tmp_H_B1           
        FLAKE1_struc(n)%flake1(t)%T_sfc           = tmp_T_sfc          
        FLAKE1_struc(n)%flake1(t)%albedo_water    = tmp_albedo_water   
        FLAKE1_struc(n)%flake1(t)%albedo_ice      = tmp_albedo_ice     
        FLAKE1_struc(n)%flake1(t)%albedo_snow     = tmp_albedo_snow    
        FLAKE1_struc(n)%flake1(t)%Q_lwa           = tmp_Q_lwa
        FLAKE1_struc(n)%flake1(t)%Q_lww           = tmp_Q_lww

        ![ 1] output variable: T_snow (unit=K). *** temperature at the air-snow interface 
        call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_LAKE_T_SNOW, &
             value = tmp_T_snow, &
             vlevel=1, unit="K", direction="-", &
             surface_type = LIS_rc%lake_index)

        ![ 2] output variable: T_ice (unit=K). *** temperature at the snow-ice interface 
        call LIS_diagnoseSurfaceOutputVar(n, t, &
             LIS_MOC_LAKE_T_ICE, value = FLAKE1_struc(n)%flake1(t)%T_ice, &
             vlevel=1, unit="K", direction="-", &
             surface_type = LIS_rc%lake_index)

        ![ 3] output variable: T_mnw (unit=K). *** mean temperature of the water column 
        call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_LAKE_T_MNW, &
             value = FLAKE1_struc(n)%flake1(t)%T_mnw, &
             vlevel=1, unit="K", direction="-", &
             surface_type = LIS_rc%lake_index)

        ![ 4] output variable: T_wML (unit=K). *** temperature of mixed layer 
        call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_LAKE_T_WML, &
             value = FLAKE1_struc(n)%flake1(t)%T_wML, &
             vlevel=1, unit="K", direction="-", &
             surface_type = LIS_rc%lake_index)

        ![ 5] output variable: T_bot (unit=K). *** temperature at the water-bottom sediment interface 
        call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_LAKE_T_BOT, &
             value = FLAKE1_struc(n)%flake1(t)%T_bot, &
             vlevel=1, unit="K", direction="-", &
             surface_type = LIS_rc%lake_index)

        ![ 6] output variable: T_b1 (unit=K). *** temperature at the bottom of the upper layer of the sediments 
        call LIS_diagnoseSurfaceOutputVar(n, t, &
             LIS_MOC_LAKE_T_B1, value = FLAKE1_struc(n)%flake1(t)%T_b1, &
             vlevel=1, unit="K", direction="-", &
             surface_type = LIS_rc%lake_index)

        ![ 7] output variable: C_T (unit=-). *** thermocline shape factor of lake
        call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_LAKE_C_T, &
             value = FLAKE1_struc(n)%flake1(t)%C_T, &
             vlevel=1, unit="-", direction="-", &
             surface_type = LIS_rc%lake_index)

        ![ 8] output variable: H_snow (unit=m). *** snow thickness above ice
        call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_SNOWDEPTH, &
             value = FLAKE1_struc(n)%flake1(t)%H_snow, &
             vlevel=1, unit="m", direction="-", &
             surface_type = LIS_rc%lake_index)

        ![ 9] output variable: H_ice (unit=m). *** ice thickness above lake
        call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_LAKE_H_ICE, &
             value = FLAKE1_struc(n)%flake1(t)%H_ice, &
             vlevel=1, unit="m", direction="-", &
             surface_type = LIS_rc%lake_index)

        ![ 10] output variable: H_ML (unit=m). *** thickness of mixed layer 
        call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_LAKE_H_ML, &
             value = FLAKE1_struc(n)%flake1(t)%H_ML, &
             vlevel=1, unit="m", direction="-", &
             surface_type = LIS_rc%lake_index)

        ![ 11] output variable: H_B1 (unit=m). *** thickness of the upper layer of bottom sediments 
        call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_LAKE_H_B1,&
             value = FLAKE1_struc(n)%flake1(t)%H_B1, &
             vlevel=1, unit="m", direction="-", &
             surface_type = LIS_rc%lake_index)

        ![ 12] output variable: T_sfc (unit=K). *** surface temperature 
        call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_AVGSURFT, &
             value = FLAKE1_struc(n)%flake1(t)%T_sfc, &
             vlevel=1, unit="K", direction="-", &
             surface_type = LIS_rc%lake_index)

        ![ 13] output variable: albedo_water (unit=-). *** water surface albedo with resect to solar radiation 
        call LIS_diagnoseSurfaceOutputVar(n, t, &
             LIS_MOC_LAKE_ALBEDO_WATER, value = &
             FLAKE1_struc(n)%flake1(t)%albedo_water, &
             vlevel=1, unit="-", direction="-", &
             surface_type = LIS_rc%lake_index)

        ![ 14] output variable: albedo_ice (unit=-). *** ice surface albedo with respect to the solar radiation 
        call LIS_diagnoseSurfaceOutputVar(n, t, &
             LIS_MOC_LAKE_ALBEDO_ICE, value = &
             FLAKE1_struc(n)%flake1(t)%albedo_ice, &
             vlevel=1, unit="-", direction="-", &
             surface_type = LIS_rc%lake_index)

        ![ 15] output variable: albedo_snow (unit=-). *** snow surface albedo with respect to the solar radiation 
        call LIS_diagnoseSurfaceOutputVar(n, t, &
             LIS_MOC_LAKE_ALBEDO_SNOW, value = &
             FLAKE1_struc(n)%flake1(t)%albedo_snow, &
             vlevel=1, unit="-", direction="-", surface_type = &
             LIS_rc%lake_index)
        ![ 16] output variable: ufr_a (unit=m/s). *** friction velocity in air 
        call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_LAKE_UFR_A, &
             value = tmp_ufr_a, &
             vlevel=1, unit="m/s", direction="-", surface_type = LIS_rc%lake_index)

        ![ 17] output variable: ufr_w (unit=m/s). *** friction velocity in surface water 
        call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_LAKE_UFR_W, &
             value = tmp_ufr_w, &
             vlevel=1, unit="m/s", direction="-", surface_type = LIS_rc%lake_index)

        ![ 18] output variable: Wconv (unit=m/s). *** convective velocity scale 
        call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_LAKE_WCONV, &
             value = tmp_Wconv, &
             vlevel=1, unit="m/s", direction="-", surface_type = LIS_rc%lake_index)

        ![ 19] output variable: Q_se (unit=W/m2). *** sensible surface heat flux 
        call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_QH, &
             value = tmp_Q_se, &
             vlevel=1, unit="W/m2", direction="UP", &
             surface_type = LIS_rc%lake_index)

        ![ 20] output variable: Q_la (unit=W/m2). *** latent surface heat flux 
        call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_QLE, &
             value = tmp_Q_la, &
             vlevel=1, unit="W/m2", direction="UP", &
             surface_type = LIS_rc%lake_index)

        ![ 21] output variable: I_w (unit=W/m2). *** radiation flux through the ice-water or air-water interface 
        call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_LAKE_I_W, &
             value = tmp_I_w, &
             vlevel=1, unit="W/m2", direction="DN", &
             surface_type = LIS_rc%lake_index)

        ![ 22] output variable: Q_lwa (unit=W/m2). *** longwave radiation flux from atmosphere to water
        call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_LWDOWNFORC, &
             value = FLAKE1_struc(n)%flake1(t)%Q_lwa, &
             vlevel=1, unit="W/m2", direction="DN", &
             surface_type = LIS_rc%lake_index)

        ![ 23] output variable: Q_lww (unit=W/m2). *** longwave radiation flux from water to atmosphere
        call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_LWUP, &
             value = FLAKE1_struc(n)%flake1(t)%Q_lww, &
             vlevel=1, unit="W/m2", direction="UP", &
             surface_type = LIS_rc%lake_index)

        ![ 24] output variable: Q_bot (unit=W/m2). *** heat flux across water-sediments boundary 
        call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_LAKE_Q_BOT, &
             value = tmp_Q_bot, &
             vlevel=1, unit="W/m2", direction="-", &
             surface_type = LIS_rc%lake_index)

        FLAKE1_struc(n)%flake1(t)%swdown = 0 
        FLAKE1_struc(n)%flake1(t)%lwdown = 0 
        FLAKE1_struc(n)%flake1(t)%wind_e = 0 
        FLAKE1_struc(n)%flake1(t)%wind_n = 0 
        FLAKE1_struc(n)%flake1(t)%tair   = 0 
        FLAKE1_struc(n)%flake1(t)%qair   = 0 
        FLAKE1_struc(n)%flake1(t)%psurf  = 0 

     enddo ! end of tile (t) loop
     FLAKE1_struc(n)%forc_count= 0 
  endif ! end of alarmcheck if 

end subroutine FLake1_main
