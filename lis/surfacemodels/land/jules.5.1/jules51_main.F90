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
! !ROUTINE: jules51_main
! \label{jules51_main}
!
! !REVISION HISTORY:
! 16 May 2016; Shugong Wang; initial implementation for JULES 4.3
! 01 Feb 2018; Shugong Wang; updated for JULES.5.1 
! 27 Feb 2020; Shugong Wang; updated for automated container 
!
! !ROUTINE: jules51_main.F90
! 
! !INTERFACE:
subroutine jules51_main(n)
   use jules_surface_mod,      only: l_aggregate
   use model_grid_mod,         only: latitude, longitude
   use model_time_mod,         only: timestep, current_time
   use jules_soil_mod,         only: dzsoil ! meter, thickness
   use jules_print_mgr,        only: jules_message, jules_print
   use LIS_coreMod,            only: LIS_rc, LIS_domain, LIS_surface
   use LIS_histDataMod
   use LIS_timeMgrMod,         only: LIS_isAlarmRinging
   use LIS_logMod,             only: LIS_logunit, LIS_endrun
   use LIS_FORC_AttributesMod
   use jules51_lsmMod
   use sf_diags_mod,           only: sf_diag
   use debug_latlon
   use c_z0h_z0m,              only: z0h_z0m
   use fluxes,                 only: emis_surft
   use csigma,                 only: sbcon
   use jules_snow_mod,         only: nsmax
   use jules_radiation_mod,    only: wght_alb
   use p_s_parms,              only: smvcst   => smvcst_soilt  , &   
                                     smvccl   => smvccl_soilt  , &   
                                     smvcwt   => smvcwt_soilt  , &   
                                     b        => bexp_soilt    , &   
                                     sathh    => sathh_soilt   , &   
                                     hcap     => hcap_soilt    , &   
                                     hcon     => hcon_soilt    , &   
                                     satcon   => satcon_soilt  , &   
                                     sthu_min => sthu_min_soilt, &
                                     albsoil  => albsoil_soilt , &
                                     catch_snow => catch_snow_surft, &
                                     infil_tile => infil_surft,     &
                                     z0_tile  => z0_surft,          &
                                     z0h_tile_bare => z0h_bare_surft
   use prognostics,            only: seed_rain,                     &
                                     smcl      => smcl_soilt       , &
                                     t_soil    => t_soil_soilt     , &
                                     lai       => lai_pft          , & 
                                     canht_ft  => canht_pft        , & 
                                     sice      => sice_surft        , &
                                     sliq      => sliq_surft       , &
                                     snowdepth => snowdepth_surft  , &
                                     tsnow     => tsnow_surft      , &
                                     ds        => ds_surft         , &
                                     rgrainl   => rgrainl_surft    , & 
                                     rho_snow_grnd => rho_snow_grnd_surft, &
                                     rho_snow  => rho_snow_surft, &
                                     canopy    => canopy_surft ,  &
                                     gc        => gc_surft,       &
                                     rgrain    => rgrain_surft,   &
                                     snow_tile => snow_surft,     &
                                     snow_grnd => snow_grnd_surft, & 
                                     tstar_tile  => tstar_surft  , &
                                     tsoil_deep => tsoil_deep_gb   !!! *** need double check 
   use ancil_info,             only: land_pts, lice_pts, soil_pts,     &
                                     lice_index, soil_index,           &
                                     l_lice_point, l_soil_point,       &
                                     land_mask,                        &
                                     tile_pts        => surft_pts    ,   &
                                     tile_index      => surft_index  , & 
                                     frac            => frac_surft   
  use forcing 
   implicit none
! !ARGUMENTS:
   integer, intent(in) :: n
!
! !DESCRIPTION:
!
!  Calls the run routines for the forcing-only option (jules51)
!
!  The arguments are:
!  \begin{description}
!  \item[n]
!   index of the nest
!  \end{description}
!EOP

   integer         :: t
   integer         :: i, k, j, m
   real            :: dt
   real            :: lat, lon
   integer         :: row, col, pft, p
   integer         :: year, month, day, hour, minute, second
   logical         :: alarmCheck
   real            :: tstar_box
   integer         :: gid, cur_grid, start_k, end_k
   logical         :: get_forcing
   real            :: lw_upward, lw_net, trad, snow_frac, albedo_land, prcp, cpcp, smc_liq, smc_ice, trad_gb_mean, smc_sat_frac

   ! EMK Additional variables
   real :: sfctmp, sfcprs, es, q2, q2sat
   real :: relsmc,smc,smcmax,smcwilt
   character*3 :: fnest

   ! check JULES alarm. If alarm is ring, run model.
   alarmCheck = LIS_isAlarmRinging(LIS_rc, "JULES.5.1 model alarm")
   if (alarmCheck) Then
      do m = 1, LIS_rc%nensem(n)
         k = 1
         do
            if ( k > LIS_rc%npatch(n,LIS_rc%lsm_index) ) then
               exit
            endif

            ! Find all tiles in current grid
            cur_grid = LIS_surface(n,LIS_rc%lsm_index)%tile(k)%index
            start_k = k
            end_k = jules51_struc(n)%jules51(start_k)%end_k 
            k = end_k + 1 
            ! reset fraction of pfts
            frac(:,:)           = 0.0
            tile_pts(:)         = 0
            tile_index(:,:)     = 0
            catch_snow(:,:)     = 0.0
            infil_tile(:,:)     = 0.0
            z0_tile(:,:)        = 0.0
            z0h_tile_bare(:,:)  = 0.0
            sice(:, :, :)       = 0.0
            sliq(:, :, :)       = 0.0
            tsnow(:, :,:)       = 0.0
            rgrainl(:, :,:)     = 0.0
            ds(:, :,:)          = 0.0
            rho_snow_grnd(:, :) = 0.0
            snowdepth(:, :)     = 0.0
            rho_snow(:, :,:)    = 0.0
            canht_ft(:, :)      = 0.0
            lai(:, :)           = 0.0
            canopy(:, :)        = 0.0
            gc(:, :)            = 0.0
            rgrain(:, :)        = 0.0
            snow_tile(:, :)     = 0.0
            snow_grnd(:, :)     = 0.0
            tstar_tile(:, :)    = 0.0
            dt = LIS_rc%nts(n)
            lat = LIS_domain(n)%grid(cur_grid)%lat
            lon = LIS_domain(n)%grid(cur_grid)%lon
            
            lat_d = lat
            lon_d = lon
            ! compute gridbox average for tstar
            tstar_box=0.0
            do t=start_k, end_k
               if ( LIS_surface(n,LIS_rc%lsm_index)%tile(t)%ensem == m ) then
                 pft = jules51_struc(n)%jules51(t)%pft 
                 jules51_struc(n)%jules51(t)%frac(pft) = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%fgrd
                 tstar_box=tstar_box+jules51_struc(n)%jules51(t)%frac(pft)*jules51_struc(n)%jules51(t)%tstar_tile(pft)
               endif
            enddo

            get_forcing = .true.
            do t=start_k, end_k
               if ( LIS_surface(n,LIS_rc%lsm_index)%tile(t)%ensem == m ) then
                  pft = jules51_struc(n)%jules51(t)%pft 
                  jules51_struc(n)%jules51(t)%frac(pft) = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%fgrd
                  jules51_struc(n)%jules51(t)%land_mask = .true.
                  jules51_struc(n)%jules51(t)%z1_tq = jules51_struc(n)%z1_tq
                  jules51_struc(n)%jules51(t)%z1_uv = jules51_struc(n)%z1_uv

                  jules51_struc(n)%jules51(t)%tstar=tstar_box

                  call tile_to_prog(n, t, pft)
                  call tile_to_ps(n, t, pft)
                  call tile_to_ancil(n, t, pft)
                  call tile_to_trifctl(n, t, pft)
                  call tile_to_bvoc_vars(n, t, pft)
                  call tile_to_jules_internal(n, t, pft)
                  call tile_to_fluxes(n, t, pft)
                  if ( get_forcing ) then
                     ! only call tile_to_forcing once
                     call tile_to_forcing(n, t)
                     call tile_to_top_pdm(n, t)
                     call tile_to_sf_diag(n, t) 
                     get_forcing = .false.
                  endif
               endif
            enddo

            ! reset time step 
            timestep = LIS_rc%tscount(n)
            current_time%year  = LIS_rc%yr
            current_time%month = LIS_rc%mo
            current_time%day   = LIS_rc%da
            current_time%time  = LIS_rc%hr * 3600 + LIS_rc%mn * 60 + LIS_rc%ss  
            
            ! model_grid_mod, ONLY: latitude, longitude
            latitude(1,1)  = LIS_domain(n)%grid(cur_grid)%lat
            longitude(1,1) = LIS_domain(n)%grid(cur_grid)%lon
            
            !------------------------------------------------------------------
            ! call jules physics 
            !------------------------------------------------------------------
            CALL lis_control() 

            do t=start_k, end_k
               if ( LIS_surface(n,LIS_rc%lsm_index)%tile(t)%ensem == m ) then
                  ! when l_aggregate == .true., pft = 1  
                  pft = jules51_struc(n)%jules51(t)%pft

                  call prog_to_tile(n, t, pft)
                  call ps_to_tile(n, t, pft)
                  call ancil_to_tile(n, t, pft)
                  call trifctl_to_tile(n,t, pft)
                  call bvoc_vars_to_tile(n, t, pft)
                  call jules_internal_to_tile(n, t, pft)
                  call fluxes_to_tile(n, t, pft)
                  call top_pdm_to_tile(n, t)
                  call sf_diag_to_tile(n, t) 

                  ![ 1] output variable: soil_temp (unit=K).
                  ! soil layer temperature
                  do i=1, jules51_struc(n)%sm_levels
                     call LIS_diagnoseSurfaceOutputVar(n, t,         &
                        LIS_MOC_SOILTEMP,                            &
                        value=jules51_struc(n)%jules51(t)%t_soil(i), &
                        vlevel=i, unit="K", direction="-",           &
                        surface_type=LIS_rc%lsm_index)
                  end do
                  ![ 2] output variable: soil_moist (unit=kg m-2).
                  ! soil layer moisture
                  do i=1, jules51_struc(n)%sm_levels
                     call LIS_diagnoseSurfaceOutputVar(n, t,       &
                        LIS_MOC_SOILMOIST,                         &
                        value=jules51_struc(n)%jules51(t)%smcl_soilt(i), &
                        vlevel=i, unit="kg m-2", direction="-",    &
                        surface_type=LIS_rc%lsm_index)
                  end do

                  ! m3/m3
                  do i=1, jules51_struc(n)%sm_levels
                     call LIS_diagnoseSurfaceOutputVar(n, t,       &
                        LIS_MOC_SOILMOIST,                         &
                        value=jules51_struc(n)%jules51(t)%smcl_soilt(i)/ &
                              (1000.0*dzsoil(i)),                  &
                        vlevel=i, unit="m^3 m-3", direction="-",   &
                        surface_type=LIS_rc%lsm_index)
                  end do
                  
                  ! m3/m3, unfrozen (liquid) soil moisture 
                  do i=1, jules51_struc(n)%sm_levels
                     ! p_s_sthu: unfrozen moisture content of each soil layer as a fraction of saturation (-)
                     ! p_s_smvcst: volumetric saturation point (m^3 m-3 of soil)
                     smc_liq = jules51_struc(n)%jules51(t)%p_s_sthu(i) * jules51_struc(n)%jules51(t)%p_s_smvcst(i) 
                     call LIS_diagnoseSurfaceOutputVar(n, t,           &
                        LIS_MOC_SMLIQFRAC,                             &
                        value=smc_liq,                                 &
                        vlevel=i, unit="m^3 m-3", direction="-",       &
                        valid_min=0.0, valid_max=1.0,                  &
                        surface_type=LIS_rc%lsm_index)
                  enddo
                  ! relative fraction, unfrozen (liquid) soil moisture 
                  do i=1, jules51_struc(n)%sm_levels
                     ! zero soil moisture is not realistic, but it could happen in JULES. 
                     smc_sat_frac = jules51_struc(n)%jules51(t)%p_s_sthu(i)+jules51_struc(n)%jules51(t)%p_s_sthf(i)
                     if (smc_sat_frac > tiny(smc_sat_frac)) then 
                        smc_liq = jules51_struc(n)%jules51(t)%p_s_sthu(i) / smc_sat_frac
                     else
                        smc_liq = 0.0
                     endif 

                     call LIS_diagnoseSurfaceOutputVar(n, t,           &
                        LIS_MOC_SMLIQFRAC,                             &
                        value=smc_liq,                                 &
                        vlevel=i, unit="-", direction="-",             &
                        valid_min=0.0, valid_max=1.0,                  &
                        surface_type=LIS_rc%lsm_index)
                  enddo
                  
                  ! m3/m3, frozen (ice) soil moisture 
                  do i=1, jules51_struc(n)%sm_levels
                     smc_ice = jules51_struc(n)%jules51(t)%p_s_sthf(i) * jules51_struc(n)%jules51(t)%p_s_smvcst(i) 
                     call LIS_diagnoseSurfaceOutputVar(n, t,           &
                        LIS_MOC_SMFROZFRAC,                            &
                        value=smc_ice,                                 &
                        vlevel=i, unit="m^3 m-3", direction="-",       &
                        valid_min=0.0, valid_max=1.0,                  &
                        surface_type=LIS_rc%lsm_index)
                  enddo

                  ! relative fraction, frozen (ice) soil moisture 
                  do i=1, jules51_struc(n)%sm_levels
                     ! zero soil moisture is not realistic, but it could happen in JULES. 
                     smc_sat_frac = jules51_struc(n)%jules51(t)%p_s_sthu(i)+jules51_struc(n)%jules51(t)%p_s_sthf(i)
                     if (smc_sat_frac > tiny(smc_sat_frac)) then
                        smc_ice = jules51_struc(n)%jules51(t)%p_s_sthf(i) / smc_sat_frac
                     else
                        smc_ice = 0.0
                     endif
                     call LIS_diagnoseSurfaceOutputVar(n, t,           &
                        LIS_MOC_SMFROZFRAC,                            &
                        value=smc_ice,                                 &
                        vlevel=i, unit="-", direction="-",             &
                        valid_min=0.0, valid_max=1.0,                  &
                        surface_type=LIS_rc%lsm_index)
                  enddo


                  do i=1, jules51_struc(n)%sm_levels
                     call LIS_diagnoseSurfaceOutputVar(n, t,           &
                        LIS_MOC_JULES_STHU,                            &
                        value=jules51_struc(n)%jules51(t)%p_s_sthu(i), &
                        vlevel=i, unit="-", direction="-",             &
                        valid_min=0.0, valid_max=1.0,                  &
                        surface_type=LIS_rc%lsm_index)
                  enddo

                  do i=1, jules51_struc(n)%sm_levels
                     call LIS_diagnoseSurfaceOutputVar(n, t,           &
                        LIS_MOC_JULES_STHF,                            &
                        value=jules51_struc(n)%jules51(t)%p_s_sthf(i), &
                        vlevel=i, unit="-", direction="-",             &
                        valid_min=0.0, valid_max=1.0,                  &
                        surface_type=LIS_rc%lsm_index)
                  enddo

                  call LIS_diagnoseSurfaceOutputVar(n, t,               &
                     LIS_MOC_JULES_STHU_MIN,                            &
                     value=jules51_struc(n)%jules51(t)%p_s_sthu_min(1), &
                     vlevel=1, unit="-", direction="-",                 &
                     valid_min=0.0, valid_max=1.0,                      &
                        surface_type=LIS_rc%lsm_index)

                  call LIS_diagnoseSurfaceOutputVar(n, t,             &
                     LIS_MOC_JULES_SMVCCL,                            &
                     value=jules51_struc(n)%jules51(t)%p_s_smvccl(1), &
                     vlevel=1, unit="m^3 m-3", direction="-",         &
                     valid_min=0.0, valid_max=1.0,                    &
                        surface_type=LIS_rc%lsm_index)

                  call LIS_diagnoseSurfaceOutputVar(n, t,             &
                     LIS_MOC_JULES_SMVCST,                            &
                     value=jules51_struc(n)%jules51(t)%p_s_smvcst(1), &
                     vlevel=1, unit="m^3 m-3", direction="-",         &
                     valid_min=0.0, valid_max=1.0,                    &
                        surface_type=LIS_rc%lsm_index)

                  call LIS_diagnoseSurfaceOutputVar(n, t,             &
                     LIS_MOC_JULES_SMVCWT,                            &
                     value=jules51_struc(n)%jules51(t)%p_s_smvcwt(1), &
                     vlevel=1, unit="m^3 m-3", direction="-",         &
                     valid_min=0.0, valid_max=1.0,                    &
                        surface_type=LIS_rc%lsm_index)

                  call LIS_diagnoseSurfaceOutputVar(n, t,  &
                     LIS_MOC_QLE,                          &
                     value=sf_diag%latent_heat(1,1),       &
                     vlevel=1,unit="W m-2",direction="UP", &
                     surface_type=LIS_rc%lsm_index)

                  call LIS_diagnoseSurfaceOutputVar(n,t,        &
                     LIS_MOC_EVAP,                              &
                     value=jules51_struc(n)%jules51(t)%fqw_1,   &
                     vlevel=1,unit="kg m-2 s-1",direction="UP", &
                     surface_type=LIS_rc%lsm_index)

                  call LIS_diagnoseSurfaceOutputVar(n, t,     &
                     LIS_MOC_QH,                              &
                     value=jules51_struc(n)%jules51(t)%ftl_1, &
                     vlevel=1,unit="W m-2",direction="UP",    &
                     surface_type=LIS_rc%lsm_index)

                  !call LIS_diagnoseSurfaceOutputVar(n, t,  &
                  !   LIS_MOC_QG,                           &
                  !   value=-gflx,                          &
                  !   vlevel=1,unit="W m-2",direction="DN", &
                  !   surface_type=LIS_rc%lsm_index)

                  call LIS_diagnoseSurfaceOutputVar(n,t,         &
                     LIS_MOC_ECANOP,                             &
                     value=jules51_struc(n)%jules51(t)%ecan,     &
                     vlevel=1,unit="kg m-2 s-1", direction="UP", &
                     surface_type=LIS_rc%lsm_index)
                  
                  call LIS_diagnoseSurfaceOutputVar(n,t,         &
                     LIS_MOC_JULES_ESOIL,                             &
                     value=jules51_struc(n)%jules51(t)%esoil,     &
                     vlevel=1,unit="kg m-2 s-1", direction="UP", &
                     surface_type=LIS_rc%lsm_index)

                  call LIS_diagnoseSurfaceOutputVar(n,t,       &
                     LIS_MOC_SUBSNOW,                          &
                     value=jules51_struc(n)%jules51(t)%ei,     &
                     vlevel=1,unit="kg m-2 s-1",direction="-", &
                     surface_type=LIS_rc%lsm_index)

                  do i=1, jules51_struc(n)%sm_levels
                     call LIS_diagnoseSurfaceOutputVar(n,t,             &
                        LIS_MOC_SOILET,                                 &
                        vlevel=i,                                       &
                        value=jules51_struc(n)%jules51(t)%ext(i),       &
                        unit="kg m-2 s-1",direction="-",                &
                        surface_type=LIS_rc%lsm_index)
                  enddo

                  call LIS_diagnoseSurfaceOutputVar(n,t,               &
                     LIS_MOC_SNOWDEPTH,                                &
                     vlevel=1,                                         &
                     value=jules51_struc(n)%jules51(t)%snowdepth(pft), &
                     unit="m", direction="-",&
                     surface_type=LIS_rc%lsm_index)

                  !(unit=mm) ***  snow water equivalent, snow_mass_ij (kg/m2) equivalent?  
                  call LIS_diagnoseSurfaceOutputVar(n, t,                &
                      LIS_MOC_SWE,                                       &
                      value =  jules51_struc(n)%jules51(t)%snow_mass_ij, &
                      vlevel=1, unit="kg m-2", direction="-",            &
                      surface_type = LIS_rc%lsm_index)
                  
                  call LIS_diagnoseSurfaceOutputVar(n, t,              &
                      LIS_MOC_WATERTABLED,                             &
                      value = JULES51_struc(n)%jules51(t)%zw,          &
                      vlevel=1, unit="m", direction="-",               &
                      surface_type = LIS_rc%lsm_index)

                  call LIS_diagnoseSurfaceOutputVar(n,t,                &
                     LIS_MOC_AVGSURFT,                                  &
                     value=jules51_struc(n)%jules51(t)%tstar_tile(pft), &
                     vlevel=1,unit="K",direction="-",                   &
                     surface_type=LIS_rc%lsm_index)

                  ! emisivity 
                  call LIS_diagnoseSurfaceOutputVar(n, t,               &
                     LIS_MOC_EMISSFORC,                                 &
                     value = jules51_struc(n)%jules51(t)%emis_tile(pft),&
                     vlevel=1, unit="-", direction="-",                 &
                     surface_type = LIS_rc%lsm_index)

                  !***  surface radiative temperature 
                  trad_gb_mean = 0.0
                  do p=1, jules51_struc(n)%jules51(t)%nsurft
                    trad_gb_mean = trad_gb_mean + frac(1,p) * jules51_struc(n)%jules51(t)%tstar_tile(p)**4.0
                  enddo
                  trad = trad_gb_mean**0.25 

                  call LIS_diagnoseSurfaceOutputVar(n, t,               &
                     LIS_MOC_RADT, &
                     value=trad, &
                     vlevel=1, unit="K", direction="-",                 &
                     surface_type=LIS_rc%lsm_index)
            
                  !***  net longwave radiation 
                  lw_upward = sbcon * jules51_struc(n)%jules51(t)%emis_tile(pft) * jules51_struc(n)%jules51(t)%tstar_tile(pft)**4 
                  lw_net    = jules51_struc(n)%jules51(t)%emis_tile(pft) * jules51_struc(n)%jules51(t)%lwdown - lw_upward 
                  call LIS_diagnoseSurfaceOutputVar(n, t,                          &
                      LIS_MOC_LWNET,vlevel=1,                                      &
                      value= lw_net,                                               &
                      unit="W m-2", direction="DN", surface_type=LIS_rc%lsm_index)

                  call LIS_diagnoseSurfaceOutputVar(n,t,        &
                     LIS_MOC_GC,                                &
                     value=jules51_struc(n)%jules51(t)%gc_surft(pft), &
                     vlevel=1,unit="m s-1",direction="-",       &
                     surface_type=LIS_rc%lsm_index)

                  call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_GS, &
                     value=jules51_struc(n)%jules51(t)%gs_gb,          &
                     vlevel=1,unit="m s-1",direction="-",           &
                     surface_type=LIS_rc%lsm_index)
                  
                  call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_JULES_STHZW, &
                     value=jules51_struc(n)%jules51(t)%sthzw,          &
                     vlevel=1,unit="-",direction="-",           &
                     surface_type=LIS_rc%lsm_index)
                  
                  call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_WATERTABLED, &
                     value=jules51_struc(n)%jules51(t)%zw,          &
                     vlevel=1,unit="m",direction="-",           &
                     surface_type=LIS_rc%lsm_index)

                  call LIS_diagnoseSurfaceOutputVar(n,t,                 &
                     LIS_MOC_ROUGHNESS,                                  &
                     value=jules51_struc(n)%jules51(t)%p_s_z0_tile(pft), &
                     vlevel=1,unit="m",direction="-",                    &
                     surface_type=LIS_rc%lsm_index)

                  call LIS_diagnoseSurfaceOutputVar(n,t,                       &
                     LIS_MOC_THERMAL_ROUGHNESS,                                &
                     value=jules51_struc(n)%jules51(t)%p_s_z0_tile(pft) * z0h_z0m(pft), &
                     vlevel=1,unit="m",direction="-",                          &
                     surface_type=LIS_rc%lsm_index)
                  
                  call LIS_diagnoseSurfaceOutputVar(n, t,                     &
                     LIS_MOC_QSB,                                             &
                     value = jules51_struc(n)%jules51(t)%sub_surf_roff,       &
                     vlevel=1, unit="kg m-2 s-1", direction="OUT", surface_type = LIS_rc%lsm_index)
                  
                  ! EMK
                  call LIS_diagnoseSurfaceOutputVar(n,t,        &
                     LIS_MOC_QS,                              &
                     value=jules51_struc(n)%jules51(t)%surf_roff*dt,   &
                     vlevel=1,unit="kg m-2",direction="OUT", &
                     surface_type=LIS_rc%lsm_index)

                  ! EMK
                  call LIS_diagnoseSurfaceOutputVar(n,t,        &
                     LIS_MOC_QSB,                              &
                     value=jules51_struc(n)%jules51(t)%sub_surf_roff*dt,   &
                     vlevel=1,unit="kg m-2",direction="OUT", &
                     surface_type=LIS_rc%lsm_index)

                  call LIS_diagnoseSurfaceOutputVar(n, t,                     &
                     LIS_MOC_QS,                                              &
                     value = jules51_struc(n)%jules51(t)%surf_roff,           &
                     vlevel=1, unit="kg m-2 s-1", direction="OUT", surface_type = LIS_rc%lsm_index)
                  
                  call LIS_diagnoseSurfaceOutputVar(n, t,                     &
                     LIS_MOC_QSM,                                             &
                     value = jules51_struc(n)%jules51(t)%snow_melt_gb,        &
                     vlevel=1, unit="kg m-2 s-1", direction="S2L", surface_type = LIS_rc%lsm_index)

                  call LIS_diagnoseSurfaceOutputVar(n, t,                     &
                     LIS_MOC_DRIP,                                            &
                     value = jules51_struc(n)%jules51(t)%tot_tfall,           &
                     vlevel=1, unit="kg m-2 s-1", direction="-", surface_type = LIS_rc%lsm_index)

                  !EMK....Fixed unit conversions.  Canopy water is already in 
                  !kg m-2.
                  call LIS_diagnoseSurfaceOutputVar(n, t,                     &
                     LIS_MOC_CANOPINT,                                        &
                     value = jules51_struc(n)%jules51(t)%canopy(pft),      &
                     vlevel=1, unit="kg m-2", direction="-", &
                     surface_type = LIS_rc%lsm_index)

                  call LIS_diagnoseSurfaceOutputVar(n, t,                     &
                     LIS_MOC_SWNET,                                           &
                     value= jules51_struc(n)%jules51(t)%sw_tile(pft),         &  
                     vlevel=1,  unit="W m-2", direction="DN", surface_type=LIS_rc%lsm_index)

                  if(pft <= jules51_struc(n)%npft) then  
                    ! vegetation surface types 
                    call LIS_diagnoseSurfaceOutputVar(n, t,                     &
                       LIS_MOC_LAI,                                             &
                       value= jules51_struc(n)%jules51(t)%lai(pft),             &  
                       vlevel=1,  unit="-", direction="-", surface_type=LIS_rc%lsm_index)
                  else 
                    ! non-vegetation surface types 
                    call LIS_diagnoseSurfaceOutputVar(n, t,                     &
                       LIS_MOC_LAI,                                             &
                       value= 0.0,             &  
                       vlevel=1,  unit="-", direction="-", surface_type=LIS_rc%lsm_index)
                  endif
                      
                  
                  call LIS_diagnoseSurfaceOutputVar(n, t,                     &
                     LIS_MOC_JULES_FSAT,                                             &
                     value= jules51_struc(n)%jules51(t)%fsat,             &  
                     vlevel=1,  unit="-", direction="-", surface_type=LIS_rc%lsm_index)
                  
                  call LIS_diagnoseSurfaceOutputVar(n, t,                     &
                     LIS_MOC_JULES_FWETL,                                             &
                     value= jules51_struc(n)%jules51(t)%fwetl,             &  
                     vlevel=1,  unit="-", direction="-", surface_type=LIS_rc%lsm_index)


! EMK BEGIN New output variables

                  ! Below code is borrowed from Noah 3.6
                  sfcprs = jules51_struc(n)%jules51(t)%psurf
                  sfctmp = jules51_struc(n)%jules51(t)%tair
                  q2 = jules51_struc(n)%jules51(t)%qair
                  es = 611.0*exp(2.501E6/461.0*(1./273.15 - 1./sfctmp))
                  q2sat = 0.622*es/(sfcprs-(1.-0.622)*es)
                  if (sfctmp .lt. &
                       jules51_struc(n)%jules51(t)%tair_agl_min) then
                     jules51_struc(n)%jules51(t)%tair_agl_min = sfctmp
                     jules51_struc(n)%jules51(t)%rhmin = q2 / q2sat
                  end if

                  call LIS_diagnoseSurfaceOutputVar(n,t,        &
                     LIS_MOC_RHMIN,                              &
                     value=jules51_struc(n)%jules51(t)%rhmin,   &
                     vlevel=1,unit="-",direction="-", &
                     surface_type=LIS_rc%lsm_index)
                  call LIS_diagnoseSurfaceOutputVar(n,t,        &
                     LIS_MOC_RHMIN,                              &
                     value=(jules51_struc(n)%jules51(t)%rhmin*100),   &
                     vlevel=1,unit="%",direction="-", &
                     surface_type=LIS_rc%lsm_index)
                  

                  ! JULES snow variables
                  if(nsmax .eq.0) then
                    ! snow ice 
                    call LIS_diagnoseSurfaceOutputVar(n, t,       &
                       LIS_MOC_SNOWICE,                         &
                       value=jules51_struc(n)%jules51(t)%sice(1,1), &
                       vlevel=1, unit="kg m-2", direction="-",    &
                       surface_type=LIS_rc%lsm_index)
                    ! snow liquid water 
                    call LIS_diagnoseSurfaceOutputVar(n, t,       &
                       LIS_MOC_SNOWLIQ,                         &
                       value=jules51_struc(n)%jules51(t)%sliq(1,1), &
                       vlevel=1, unit="kg m-2", direction="-",    &
                       surface_type=LIS_rc%lsm_index)
                    ! snow temperature
                    call LIS_diagnoseSurfaceOutputVar(n, t,       &
                       LIS_MOC_SNOWT,                         &
                       value=jules51_struc(n)%jules51(t)%tsnow(1,1), &
                       vlevel=1, unit="K", direction="-",    &
                       surface_type=LIS_rc%lsm_index)
                    ! bulk snow density 
                    call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_SNOWDENSITY, &
                        value =  jules51_struc(n)%jules51(t)%rho_snow_grnd(pft),&
                        vlevel=1, unit="kg m-3", direction="-", surface_type = LIS_rc%lsm_index)
                    
                    ! bulk snow grain size 
                    call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_SNOWGRAIN, &
                        value =  jules51_struc(n)%jules51(t)%rgrain(pft),&
                        vlevel=1, unit="micron", direction="-", surface_type = LIS_rc%lsm_index)


                    call LIS_diagnoseSurfaceOutputVar(n,t,               &
                       LIS_MOC_SNOWTHICK,                                &
                       vlevel=1,                                         &
                       value=jules51_struc(n)%jules51(t)%snowdepth(pft), &
                       unit="m", direction="-",&
                       surface_type=LIS_rc%lsm_index)

                  else
                    do i=1, jules51_struc(n)%nsmax
                      ! snow ice 
                       call LIS_diagnoseSurfaceOutputVar(n, t,       &
                          LIS_MOC_SNOWICE,                         &
                          value=jules51_struc(n)%jules51(t)%sice(pft,i), &
                          vlevel=i, unit="kg m-2", direction="-",    &
                          surface_type=LIS_rc%lsm_index)
                      ! snow liquid water 
                       call LIS_diagnoseSurfaceOutputVar(n, t,       &
                          LIS_MOC_SNOWLIQ,                         &
                          value=jules51_struc(n)%jules51(t)%sliq(pft,i), &
                          vlevel=i, unit="kg m-2", direction="-",    &
                          surface_type=LIS_rc%lsm_index)
                      ! snow temperature K 
                       call LIS_diagnoseSurfaceOutputVar(n, t,            &
                          LIS_MOC_SNOWTPROF,                              &
                          value=jules51_struc(n)%jules51(t)%tsnow(pft,i), &
                          vlevel=i, unit="K", direction="-",              &
                          surface_type=LIS_rc%lsm_index)
                      ! snow density 
                      call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_SNOWDENSITY, &
                        value = jules51_struc(n)%jules51(t)%rho_snow(pft,i),&
                        vlevel=i, unit="kg m-3", direction="-", surface_type = LIS_rc%lsm_index)
                      
                      ! snow grain size 
                      call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_SNOWGRAIN, &
                        value =  jules51_struc(n)%jules51(t)%rgrainl(pft,i),&
                        vlevel=i, unit="micron", direction="-", surface_type = LIS_rc%lsm_index)

                      ! thickness of snow layers  
                      call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_SNOWTHICK, &
                        value = jules51_struc(n)%jules51(t)%ds(pft,i),& 
                        vlevel=i, unit="m", direction="-", surface_type = LIS_rc%lsm_index)
                    end do
                  endif

                  !snow frac according to the algorithm in JULES output interface 
                  snow_frac = 0.0
                  if(l_aggregate) then
                    if(jules51_struc(n)%jules51(t)%snow_tile(1)+jules51_struc(n)%jules51(t)%snow_grnd(1)>1.0) then
                      snow_frac = 1.0
                    endif
                  else
                    if(jules51_struc(n)%jules51(t)%snow_tile(pft)+jules51_struc(n)%jules51(t)%snow_grnd(pft)>1.0) then
                      snow_frac = jules51_struc(n)%jules51(t)%frac(pft)
                    endif
                  endif
                  call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_SNOWCOVER,             &
                    value = snow_frac,                                                   &
                    vlevel=1, unit="-", direction="-", surface_type = LIS_rc%lsm_index)
                  call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_SNOWCOVER,             &
                    value = (snow_frac*100.0),                                           &
                    vlevel=1, unit="%", direction="-", surface_type = LIS_rc%lsm_index)

                  prcp = jules51_struc(n)%jules51(t)%rainf/jules51_struc(n)%forc_count
                  
                  if (jules51_struc(n)%jules51(t)%tair  .lt. 273.16) then 
                     call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_SNOWF,value=prcp,       &    
                          vlevel=1,unit="kg m-2 s-1",direction="DN",surface_type=LIS_rc%lsm_index)
                     call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_RAINF,value=0.0,        &    
                          vlevel=1,unit="kg m-2 s-1",direction="DN",surface_type=LIS_rc%lsm_index)
                     call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_SNOWF,value=prcp*dt,    &    
                          vlevel=1,unit="kg m-2",direction="DN",surface_type=LIS_rc%lsm_index)
                     call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_RAINF,value=0.0,        &    
                          vlevel=1,unit="kg m-2",direction="DN",surface_type=LIS_rc%lsm_index)
                  else 
                     call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_SNOWF,value=0.0,        &    
                          vlevel=1,unit="kg m-2 s-1",direction="DN",surface_type=LIS_rc%lsm_index)
                     call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_RAINF,value=prcp,       &    
                          vlevel=1,unit="kg m-2 s-1",direction="DN",surface_type=LIS_rc%lsm_index)
                     call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_SNOWF,value=0.0,        &    
                          vlevel=1,unit="kg m-2",direction="DN",surface_type=LIS_rc%lsm_index)
                     call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_RAINF,value=prcp*dt,    &    
                          vlevel=1,unit="kg m-2",direction="DN",surface_type=LIS_rc%lsm_index)
                  endif

                  cpcp   = jules51_struc(n)%jules51(t)%rainf_c/jules51_struc(n)%forc_count

                  call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_CRAINFFORC,value=cpcp,      &
                       vlevel=1,unit="kg m-2 s-1",direction="DN",&
                       surface_type=LIS_rc%lsm_index,valid_min=0.0,valid_max=0.02)
                  call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_CRAINFFORC,value=cpcp*dt,   &
                       vlevel=1,unit="kg m-2",direction="DN",surface_type=LIS_rc%lsm_index)
                  
                  call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_LSRAINF,value=ls_rain_ij(1,1),      &
                       vlevel=1,unit="kg m-2 s-1",direction="DN",&
                       surface_type=LIS_rc%lsm_index,valid_min=0.0,valid_max=0.02)
                  call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_CRAINF,value=con_rain_ij(1,1),      &
                       vlevel=1,unit="kg m-2 s-1",direction="DN",&
                       surface_type=LIS_rc%lsm_index,valid_min=0.0,valid_max=0.02)
                  call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_LSSNOWF,value=ls_snow_ij(1,1),      &
                       vlevel=1,unit="kg m-2 s-1",direction="DN",&
                       surface_type=LIS_rc%lsm_index,valid_min=0.0,valid_max=0.02)
                  call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_CSNOWF,value=con_snow_ij(1,1),      &
                       vlevel=1,unit="kg m-2 s-1",direction="DN",&
                       surface_type=LIS_rc%lsm_index,valid_min=0.0,valid_max=0.02)
                   
                  call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_LSRAINF,value=ls_rain_ij(1,1) * dt,      &
                       vlevel=1,unit="kg m-2",direction="DN",&
                       surface_type=LIS_rc%lsm_index,valid_min=0.0,valid_max=200.0)
                  call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_CRAINF,value=con_rain_ij(1,1) * dt,      &
                       vlevel=1,unit="kg m-2",direction="DN",&
                       surface_type=LIS_rc%lsm_index,valid_min=0.0,valid_max=200.0)
                  call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_LSSNOWF,value=ls_snow_ij(1,1) * dt,      &
                       vlevel=1,unit="kg m-2",direction="DN",&
                       surface_type=LIS_rc%lsm_index,valid_min=0.0,valid_max=200.0)
                  call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_CSNOWF,value=con_snow_ij(1,1) * dt,      &
                       vlevel=1,unit="kg m-2",direction="DN",&
                       surface_type=LIS_rc%lsm_index,valid_min=0.0,valid_max=200.0)

                  ! albedo land based on JULES algorithm  Shugong Wang 05/02/2018
                  ! Calculate the albedo as used in subroutine control when calculating the net
                  ! shortwave on tiles. Here we take the average of diffuse albedos in VIS and NIR
                  albedo_land = wght_alb(1) * jules51_struc(n)%jules51(t)%land_albedo(1) &
                              + wght_alb(2) * jules51_struc(n)%jules51(t)%land_albedo(2) &
                              + wght_alb(3) * jules51_struc(n)%jules51(t)%land_albedo(3) &
                              + wght_alb(4) * jules51_struc(n)%jules51(t)%land_albedo(4) 

                  call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_ALBEDO,                &
                    value = albedo_land,                                                 &
                    vlevel=1, unit="-", direction="-", surface_type = LIS_rc%lsm_index)
                  call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_ALBEDO,                &
                    value = albedo_land*100.0,                                           &
                    vlevel=1, unit="%", direction="-", surface_type = LIS_rc%lsm_index)
                  
                  call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_ALBEDO_DIR_V,            &
                    value = jules51_struc(n)%jules51(t)%land_albedo(1),                    &
                    vlevel=1, unit="-", direction="-", surface_type = LIS_rc%lsm_index)
                  call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_ALBEDO_DIF_V,            &
                    value = jules51_struc(n)%jules51(t)%land_albedo(2),                    &
                    vlevel=1, unit="-", direction="-", surface_type = LIS_rc%lsm_index)
                  call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_ALBEDO_DIR_N,            &
                    value = jules51_struc(n)%jules51(t)%land_albedo(3),                    &
                    vlevel=1, unit="-", direction="-", surface_type = LIS_rc%lsm_index)
                  call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_ALBEDO_DIF_N,            &
                    value = jules51_struc(n)%jules51(t)%land_albedo(4),                    &
                    vlevel=1, unit="-", direction="-", surface_type = LIS_rc%lsm_index)
                  
                  call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_COSZENFORC,            &
                    value = jules51_struc(n)%jules51(t)%p_s_cosz,                        &
                    vlevel=1, unit="-", direction="-", surface_type = LIS_rc%lsm_index)

!                 ! JULES GPP and NPP
                  ! output variable: gpp (unit=g m-2 s-1 ). ***  net instantaneous assimilation of carbon, JULES unit: (kg C/m2/s) 
                  ! (g m-2 s-1) = 1000.0 * (kg m-2 s-1)
                  call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_GPP, value = 1000.0 * jules51_struc(n)%jules51(t)%gpp, &
                                              vlevel=1, unit="g m-2 s-1", direction="IN", surface_type = LIS_rc%lsm_index)
            
                  ! output variable: npp (unit=g m-2 s-1). ***  net primary productivity of carbon, JULES unit: (kg C/m2/s)
                  ! (g m-2 s-1) = 1000.0 * (kg m-2 s-1)
                  call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_NPP, value = 1000.0 * jules51_struc(n)%jules51(t)%npp, &
                                              vlevel=1, unit="g m-2 s-1", direction="IN", surface_type = LIS_rc%lsm_index)
!                ! End JULES GPP and NPP 

!                   !EMK Below code is borrowed from Noah 3.6
                  do i=1, jules51_struc(n)%sm_levels
                     smc = jules51_struc(n)%jules51(t)%smcl_soilt(i)/ &
                          (1000.0*dzsoil(i))
                     smcwilt = jules51_struc(n)%jules51(t)%sm_wilt
                     smcmax = jules51_struc(n)%jules51(t)%sm_sat
                     if ((smcmax - smcwilt) > 0) then
                        relsmc = (smc - smcwilt)/(smcmax - smcwilt)
                        if (relsmc > 1.0) then
                           relsmc = 1.0
                        else if (relsmc < 0.01) then
                           relsmc = 0.01
                        end if
                     else
                        relsmc = LIS_rc%udef
                     end if
                      call LIS_diagnoseSurfaceOutputVar(n,t,        &
                           LIS_MOC_RELSMC,                              &
                           value=relsmc,&
                           vlevel=i,unit="-",direction="-", &
                           surface_type=LIS_rc%lsm_index)
                      if (relsmc .ne. LIS_rc%udef) then
                         relsmc = relsmc*100
                      end if
                      call LIS_diagnoseSurfaceOutputVar(n,t,        &
                           LIS_MOC_RELSMC,                              &
                           value=relsmc,&
                           vlevel=i,unit="%",direction="-", &
                           surface_type=LIS_rc%lsm_index)
                   end do ! i
                  
! EMK END New output variables

                  ! reset forcing variables to 0 for accumulation
                  jules51_struc(n)%jules51(t)%tair    = 0.0
                  jules51_struc(n)%jules51(t)%psurf   = 0.0
                  jules51_struc(n)%jules51(t)%wind_e  = 0.0
                  jules51_struc(n)%jules51(t)%wind_n  = 0.0
                  jules51_struc(n)%jules51(t)%qair    = 0.0
                  jules51_struc(n)%jules51(t)%swdown  = 0.0
                  jules51_struc(n)%jules51(t)%lwdown  = 0.0
                  jules51_struc(n)%jules51(t)%rainf   = 0.0
                  jules51_struc(n)%jules51(t)%rainf_c = 0.0
                  jules51_struc(n)%jules51(t)%snowf   = 0.0
               endif
            enddo
         enddo
      enddo
      jules51_struc(n)%forc_count = 0
   endif

   ! EMK...See if jules51_struc(n)%jules51(t)%tair_agl_min needs to be reset
   ! for calculating RHMin.  
   write(fnest,'(i3.3)') n
   alarmCheck = LIS_isAlarmRinging(LIS_rc, &
        "JULES.5.1 RHMin alarm "//trim(fnest))
   if (alarmCheck) then
      write(LIS_logunit,*)'[INFO] Resetting tair_agl_min for RHMin calculation'
     do t = 1,LIS_rc%npatch(n,LIS_rc%lsm_index)
        jules51_struc(n)%jules51(t)%tair_agl_min = 999.
     end do
  end if

end subroutine jules51_main
