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
! !ROUTINE: snowmodel_main
! \label{snowmodel_main}
! 
! !REVISION HISTORY:
!   14 Apr 2020: Kristi Arsenault; Add G. Liston's SnowModel 
!
! !INTERFACE:
subroutine snowmodel_main(n)
! !USES:
  use LIS_coreMod
  use LIS_timeMgrMod,    only : LIS_isAlarmRinging
  use LIS_logMod,        only : LIS_logunit, LIS_endrun
  use LIS_histDataMod
  use LIS_FORC_AttributesMod
  use LIS_mpiMod
  use LIS_metforcingMod, only : LIS_forc

  use snowmodel_module
  use snowmodel_lsmMod
  use snowmodel_inc
  use snowmodel_vars
  use micromet_module

!
! !DESCRIPTION:
!  This is the entry point for calling the SnowModel physics.
!  This routine calls the {\tt ENBAL} routine that performs the
!  land surface computations, to solve for water and energy equations.
!  For documentation of the {\tt ENBAL} and routines from the
!   SnowModel, please see:
!   ftp://gliston.cira.colostate.edu/SnowModel/code/ 
! 
!  The arguments are: 
!  \begin{description}
!  \item[n]
!   index of the nest
!  \end{description}
!EOP
  implicit none

  integer, intent(in)   :: n

  integer           :: t, ierr
  integer           :: row, col, c, r
  double precision  :: xmn_part  ! center x of local LL starting point
  double precision  :: ymn_part  ! center y of local LL starting point

  integer           :: J_day     ! model Julian day, actually day-of-year
  integer           :: nveg
  real              :: rad2deg
  real, allocatable :: sphm_grid(:,:)

  real              :: delta_topo
  real              :: T_lapse_rate, Td_lapse_rate
  real              :: precip_lapse_rate
  real              :: forcelev
  real              :: forest_frac, trans_veg

  ! Vegetation impact on winds:
  real              :: beta, a
  real              :: veg_ht, canopy_windwt

! Define the Stefan Boltzmann constant.
  real, parameter :: Stef_Boltz = 5.6696e-8

  character*3       :: fnest
  logical           :: isTimeToRunCheck

! --------------------------------------------

  write(fnest,'(i3.3)') n
  isTimeToRunCheck = LIS_isAlarmRinging(LIS_rc, "SnowModel model alarm "//trim(fnest) ) 

  if(isTimeToRunCheck .neqv. .true.) return

    write(LIS_logunit,*) '[INFO] Call to the SnowModel Main routine ...'

    ! Using local parallel subdomain starting index values:
    ! LIS_ews_halo_ind(n,LIS_localPet+1) -- defined in LIS_coreMod.F90 ...
    xmn_part = xmn + deltax * ( real(LIS_ews_halo_ind(n,LIS_localPet+1)) - 1.0 )
    ymn_part = ymn + deltay * ( real(LIS_nss_halo_ind(n,LIS_localPet+1)) - 1.0 )

    iter_start = 1
    snowmodel_struc(n)%iter = snowmodel_struc(n)%iter + 1
    iter = snowmodel_struc(n)%iter

    write(LIS_logunit,*) "[INFO] Print Snowmodel iter counter: ",iter

    ! Determine source of calls to MicroMet-based routines:
    if( snowmodel_struc(n)%sm_micromet_opt == "LIS" ) then

      ! Calculate what the current simulation date should be.
      call get_model_time(iyear_init,imonth_init,iday_init,&
                          xhour_init,iter,dt,iyear,imonth,iday,xhour,J_day)

        if( LIS_masterproc ) then
           write(LIS_logunit,151) iyear,imonth,iday,xhour
        endif
  151   format('                   [INFO] SNOWMODEL TIME =',&
               i5,2i4,f6.1)

      ! _________________________________________________________

      ! Convert local LIS npatch array to the the SnowModel nx,ny grids
      if( LIS_masterproc ) then
        write(LIS_logunit,*) "[INFO] Replacing MicroMet grid forcings with LIS forcing "
        write(LIS_logunit,*) "  main:", LIS_rc%mo, LIS_rc%da, LIS_rc%hr, snowmodel_struc(n)%forc_count
      endif

      ! Constants required for wind direction and RH calculations:
      rad2deg = 180.0 / (2.0 * acos(0.0))  ! pi = (2.0 * acos(0.0))
      windspd_flag = 0.0

      allocate( sphm_grid(LIS_rc%lnc(n),LIS_rc%lnr(n)) )

      do t = 1,LIS_rc%npatch(n,LIS_rc%lsm_index)
         col = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%col
         row = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%row

         ! Translate LIS-based metforcing into SnowModel required inputs:
         tair_grid(col,row) = snowmodel_struc(n)%sm(t)%tair / snowmodel_struc(n)%forc_count
         Qsi_grid(col,row) = snowmodel_struc(n)%sm(t)%swdown / snowmodel_struc(n)%forc_count
         Qli_grid(col,row) = snowmodel_struc(n)%sm(t)%lwdown / snowmodel_struc(n)%forc_count
         sfc_pressure(col,row) = snowmodel_struc(n)%sm(t)%psurf / snowmodel_struc(n)%forc_count

         ! Total precip (in meters)
         prec_grid(col,row) = snowmodel_struc(n)%sm(t)%rainf / snowmodel_struc(n)%forc_count

         ! Convert the precipitation values from mm to m swe.  Also, make
         !   sure the interpolation has not created any negetive
         !   precipitation values.
         prec_grid(col,row) = (prec_grid(col,row) / 1000.0) * 3600 * snowmodel_struc(n)%forc_count
         prec_grid(col,row) = max(0.0,prec_grid(col,row))

         ! U- and V-gridded winds:
         uwind_grid(col,row) = snowmodel_struc(n)%sm(t)%uwind &
                             / snowmodel_struc(n)%forc_count
         vwind_grid(col,row) = snowmodel_struc(n)%sm(t)%vwind &
                             / snowmodel_struc(n)%forc_count

         ! Wind speed grid:
         windspd_grid(col,row) = sqrt(uwind_grid(col,row)**2 + vwind_grid(col,row)**2)
         ! Max windspeed flag:
         windspd_flag = max(windspd_flag,windspd_grid(col,row)) ! Handled below in "micromet_wind" call

         ! Wind direction:
         ! Some compilers do not allow both u and v to be 0.0 in
         !  the atan2 computation.
         if( abs(uwind_grid(col,row)).lt.1e-10 ) then 
           uwind_grid(col,row) = 1e-10
         endif 
         winddir_grid(col,row) = rad2deg * atan2(uwind_grid(col,row),vwind_grid(col,row))
         if( winddir_grid(col,row).ge.180.0 ) then
           winddir_grid(col,row) = winddir_grid(col,row) - 180.0
         else
           winddir_grid(col,row) = winddir_grid(col,row) + 180.0
         endif

         ! Specific humidity:
         sphm_grid(col,row) = snowmodel_struc(n)%sm(t)%qair / snowmodel_struc(n)%forc_count

      end do


      ! Call functions to derive additional forcing inputs to SnowModel:
      call micromet_relhum( LIS_rc%lnc(n),LIS_rc%lnr(n), & 
                            tair_grid, sphm_grid, sfc_pressure, rh_grid )

      call micromet_snowfall(LIS_rc%lnc(n),LIS_rc%lnr(n), &
                             snowfall_frac, prec_grid, tair_grid, sprec)

      deallocate( sphm_grid )

      ! Calculate the forest lai for each of the five forest types, and
      !   for this day of the simulation (the lai varies seasonally for
      !   the case of deciduous trees).
      call get_lai(J_day,forest_LAI)


      ! Account for vegetation impacts on SWdown, LWdown, and winds :
      do t = 1,LIS_rc%npatch(n,LIS_rc%lsm_index)
         col = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%col
         row = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%row

         if( vegtype(col,row).le.5.0 ) then
           if( calc_subcanopy_met.eq.1.0 ) then

             nveg = nint(vegtype(col,row))

             ! [SWdown] Incoming solar radiation veg-modifications:
             ! Define the forest-canopy transmissivity.  0.71 provided a
             !   best-fit to the observations, when averaged over the two years
             !   of hourly data.
             trans_veg = exp((- 0.71) * forest_LAI(nveg))

             ! Account for any gaps in the forest canopy that will allow
             !   direct incoming solar radiation to reach the snow surface.
             trans_veg = gap_frac * (1.0 - trans_veg) + trans_veg
             Qsi_grid(col,row) = trans_veg * Qsi_grid(col,row)

             ! [LWdown] Modify the incoming longwave radiation for the forest canopy.
             if( forest_LAI(nveg).lt.0.2 ) then
                forest_frac = 0.5 * forest_LAI(nveg)
             else
                forest_frac = &
                    min(1.0,0.55 + 0.29 * log(forest_LAI(nveg)))
             endif
             Qli_grid(col,row) = Qli_grid(col,row) * (1.0 - forest_frac) + &
                  (Stef_Boltz * Tair_grid(col,row)**4) * forest_frac

             !i [Wind] Further modify the wind speed to account for forest canopies.
             ! Define the canopy wind-weighting factor.  Assume z=0.6*canopy_ht,
             !   and the canopy_ht equals the vegetation snow-holding depth.
              beta = 0.9
              veg_ht = vegsnowd_xy(col,row)
              a = beta * forest_LAI(nveg)
              canopy_windwt = exp((- a)*(1.0 - (0.6*veg_ht)/veg_ht))
              windspd_grid(col,row) = canopy_windwt * windspd_grid(col,row)

           endif
         endif
      enddo

      ! Need to deterime max windspeed for entire domain:
#if (defined SPMD)
      call MPI_Barrier(LIS_MPI_COMM, ierr)
      call MPI_ALLREDUCE(windspd_flag, snowmodel_struc(n)%windspdflg_glb, 1,&
           MPI_REAL, MPI_MAX,&
           LIS_mpi_comm, ierr)
      windspd_flag = snowmodel_struc(n)%windspdflg_glb
#endif

#if 0
      ! Make the topographic calculations required by the wind and solar
      !   radiation models.  These calculations are not fixed in time
      !   because as the snow depth evolves it modifies the "topography".
      call topo_data(LIS_rc%lnc(n),LIS_rc%lnr(n),deltax,deltay,topo,&
                     curvature,terrain_slope,slope_az,curve_len_scale)

      ! Calculate the temperature and dew-point lapse rates to be used in
      !   the interpolations.
      call get_lapse_rates(imonth,iday,T_lapse_rate,&
               Td_lapse_rate,xlat_grid(1,1),lapse_rate_user_flag,&
               precip_lapse_rate,iprecip_lapse_rate_user_flag)
#endif

#if 0
      ! WIND SPEED AND DIRECTION.
      if (i_wind_flag.eq.1) then
        call micromet_wind( LIS_rc%lnc(n),LIS_rc%lnr(n),deltax,deltay, &
                         uwind_grid,vwind_grid,slopewt,curvewt,&
                         curvature,slope_az,terrain_slope,windspd_grid,&
                         winddir_grid,windspd_flag,winddir_flag,windspd_min,&
                         vegtype,forest_LAI,calc_subcanopy_met,vegsnowd_xy,&
                         topo,wind_lapse_rate,curve_len_scale)
      else
      ! Need to deterime max windspeed for entire domain:
#if (defined SPMD)
      call MPI_Barrier(LIS_MPI_COMM, ierr)
      call MPI_ALLREDUCE(windspd_flag, snowmodel_struc(n)%windspdflg_glb, 1,&
           MPI_REAL, MPI_MAX,&
           LIS_mpi_comm, ierr)
      windspd_flag = snowmodel_struc(n)%windspdflg_glb
#endif
      endif
#endif

#if 0
      ! SOLAR RADIATION.
      if (i_solar_flag.eq.1) then
        call solar( LIS_rc%lnc(n),LIS_rc%lnr(n),&
             xhour,J_day,topo,rh_grid,Tair_grid,&
             xlat_grid,Qsi_grid,slope_az,terrain_slope,dt,vegtype,&
             forest_LAI,T_lapse_rate,Td_lapse_rate,&
             calc_subcanopy_met,gap_frac,cloud_frac_factor,UTC_flag,&
             xlon_grid,cloud_frac_grid)
      endif
#endif 

! -------------------------------------------------------------------  

! *** Read in SnowModel preprocessed forcing, which can be directly
!      fed into SnowModel's MicroMet routines:

   elseif( snowmodel_struc(n)%sm_micromet_opt == "SnowModel" .and. &
           run_micromet.eq.1.0 ) then

      ! SnowModel calls: Distribute the meteorological station data.
      CALL MICROMET_CODE(LIS_rc%lnc(n),LIS_rc%lnr(n),xmn_part,ymn_part,deltax,deltay, &
         iyear_init,imonth_init,iday_init,xhour_init,dt,undef,&
         ifill,iobsint,dn,iter,curve_len_scale,slopewt,curvewt,&
         topo,curvature,terrain_slope,slope_az,Tair_grid,&
         rh_grid,uwind_grid,vwind_grid,Qsi_grid,prec_grid,&
         i_tair_flag,i_rh_flag,i_wind_flag,i_solar_flag,&
         i_prec_flag,isingle_stn_flag,igrads_metfile,&
         windspd_grid,winddir_grid,windspd_flag,winddir_flag,&
         sprec,windspd_min,Qli_grid,i_longwave_flag,vegtype,&
         forest_LAI,iyear,imonth,iday,xhour,corr_factor,&
         icorr_factor_index,lapse_rate_user_flag,&
         iprecip_lapse_rate_user_flag,use_shortwave_obs,&
         use_longwave_obs,use_sfc_pressure_obs,sfc_pressure,&
         run_enbal,run_snowpack,calc_subcanopy_met,vegsnowd_xy,&
         gap_frac,cloud_frac_factor,barnes_lg_domain,n_stns_used,&
         k_stn,xlat_grid,xlon_grid,UTC_flag,icorr_factor_loop,&
         snowmodel_line_flag,xg_line,yg_line,irun_data_assim,&
         wind_lapse_rate,iprecip_scheme,cf_precip_flag,cf_precip,&
         cloud_frac_grid,snowfall_frac,seaice_run,metforce_opt)

   else
      write(LIS_logunit,*) "[ERR] Incorrect option set for "
      write(LIS_logunit,*) "[ERR]  SnowModel MicroMet input source: "
      write(LIS_logunit,*) "[ERR]  See documentation in configs/lis.config.adoc "
      write(LIS_logunit,*) "[ERR]  for option details. "
      call LIS_endrun
   endif

   ! Perform a surface energy balance over the domain.
   if (run_enbal.eq.1.0) then
       CALL ENBAL_CODE(LIS_rc%lnc(n),LIS_rc%lnr(n),Tair_grid,uwind_grid,sfc_pressure,&
         vwind_grid,rh_grid,Tsfc,Qsi_grid,Qli_grid,Qle,Qh,Qe,&
         Qc,Qm,e_balance,Qf,snow_d,ht_windobs,icond_flag,&
         albedo,snow_z0,veg_z0,vegtype,undef,albedo_snow_forest,&
         albedo_snow_clearing,albedo_glacier,snod_layer,T_old,&
         gamma,KK)
   endif

   ! Evolve the snowpack according to the defined melt and
   !   precipitation inputs.
   if (run_snowpack.eq.1.0) then
      CALL SNOWPACK_CODE(LIS_rc%lnc(n),LIS_rc%lnr(n),Tair_grid,rh_grid,ro_nsnow,&
         dt,swe_depth,Tsfc,snow_d,prec_grid,runoff,Qm,rain,&
         sprec,iter,w_balance,sum_prec,sum_runoff,xro_snow,&
         undef,ro_snow,ro_snow_grid,soft_snow_d,sum_sprec,&
         snow_depth,windspd_grid,Qsi_grid,sum_Qcs,canopy_int,&
         Qcs,vegtype,forest_LAI,albedo,glacier_melt,&
         canopy_unload,sum_unload,sum_glacmelt,run_snowtran,&
         swemelt,d_canopy_int,sum_d_canopy_int,snow_d_init,&
         sfc_pressure,Qe,sfc_sublim_flag,sum_sfcsublim,&
         sum_swemelt,corr_factor,icorr_factor_index,swesublim,&
         swe_depth_old,canopy_int_old,KK,max_layers,melt_flag,&
         ro_snowmax,tsls_threshold,dz_snow_min,tslsnowfall,&
         change_layer,snod_layer,swed_layer,ro_layer,T_old,gamma,&
         multilayer_snowpack,seaice_run,seaice_conc,ht_windobs,&
         windspd_2m_grid,diam_layer,flux_layer)
   endif

   ! Run the blowing-snow model.
   if (run_snowtran.eq.1.0) then
      CALL SNOWTRAN_CODE(bc_flag,bs_flag,C_z,&
         conc_salt,deltax,deltay,dh_salt,dh_salt_u,dh_salt_v,&
         dh_susp,dh_susp_u,dh_susp_v,dt,dz_susp,fall_vel,fetch,&
         gravity,h_const,h_star,ht_rhobs,ht_windobs,index_ue,&
         index_uw,index_vn,index_vs,iter,LIS_rc%lnc(n),LIS_rc%lnr(n),pi,Qsalt,Qsalt_max,&
         Qsalt_maxu,Qsalt_maxv,Qsalt_u,Qsalt_v,Qsubl,Qsusp,&
         Qsusp_u,Qsusp_v,rh_grid,ro_air,ro_snow,ro_water,snow_d,&
         snow_d_init,snow_z0,soft_snow_d,sprec,sum_glacmelt,&
         subgrid_flag,wbal_salt,wbal_susp,wbal_qsubl,sum_sprec,&
         tabler_ee,tabler_ne,tabler_nn,tabler_nw,tabler_se,&
         tabler_ss,tabler_sw,tabler_ww,tair_grid,topo,topo_land,&
         topoflag,twolayer_flag,Up_const,Ur_const,Utau,&
         Utau_t,uwind_grid,veg_z0,vegsnowd_xy,vegtype,vonKarman,&
         vwind_grid,wind_min,winddir_flag,winddir_grid,&
         windspd_flag,windspd_grid,xmu,z_0,ztop_susp,max_iter,&
         run_enbal,run_snowpack,wbal_subgrid,sum_qsubl,sum_trans,&
         swe_depth,snow_depth,ro_snow_grid,sum_prec,sum_runoff,&
         sum_Qcs,canopy_int,w_balance,sum_sfcsublim,tabler_dir,&
         slope_adjust,Utau_t_const,Utau_t_flag,ro_soft_snow_old,&
         ro_soft_snow,ro_nsnow,prec_grid,Qcs,runoff,d_canopy_int,&
         glacier_melt,swe_depth_old,swesublim,canopy_unload,&
         canopy_int_old,iter_start,multilayer_snowpack,swed_layer,&
         KK,snod_layer,ro_layer,curve_lg_scale_flag,curve_wt_lg,&
         seaice_run,seaice_conc,tslsnowfall,T_old,tsls_threshold,&
         curve_len_scale,Tabler_1_flag,Tabler_2_flag,undef,&
         tabler_sfc_path_name,output_path_wo_assim,&
         output_path_wi_assim,icorr_factor_loop,windspd_2m_grid,&
         Qsubl_depth)
   endif


#if 0
   ! Save the outputs from the SNOWPACK and SNOWTRAN routines
   if (print_user.eq.1.0) then
      CALL OUTPUTS_USER(LIS_rc%lnc(n),LIS_rc%lnr(n),iter,Tair_grid,rh_grid,&
         uwind_grid,vwind_grid,windspd_grid,winddir_grid,&
         Qsi_grid,Qli_grid,prec_grid,Tsfc,Qle,Qh,Qe,Qc,Qm,Qf,&
         e_balance,snow_depth,xro_snow,swe_depth,ro_nsnow,&
         runoff,rain,sprec,sum_prec,sum_runoff,w_balance,&
         snow_d,topo_land,wbal_qsubl,sum_sprec,wbal_salt,&
         wbal_susp,ro_snow_grid,sum_Qcs,canopy_int,Qcs,&
         iyear,imonth,iday,xhour,undef,deltax,xmn,ymn,&
         wbal_subgrid,canopy_unload,sum_qsubl,sum_trans,&
         sum_unload,sum_glacmelt,glacier_melt,swemelt,&
         sfc_pressure,sum_swemelt,albedo,nrecs_max,&
         icorr_factor_loop,swesublim,vegtype,iter_start,&
         seaice_run,print_inc,cloud_frac_grid,&
         output_path_wo_assim,output_path_wi_assim,print_var,&
         print_outvars,Qsubl_depth)
   endif
#endif

!- Write out final 1D variables to LIS output files:
   do t = 1,LIS_rc%npatch(n,LIS_rc%lsm_index)

      ! Need to convert the SnowModel nx,ny grids to 
      !  local LIS npatch array ...
      col = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%col
      row = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%row

      !-- Meteorological forcing fields --

      ! Air temperature (units: K; Tair_grid)
      !  Snowmodel units: degrees C (tair);
      snowmodel_struc(n)%sm(t)%tair = &
            tair_grid(col,row)

      ! Downward radiation components:
      snowmodel_struc(n)%sm(t)%swdown = &
            Qsi_grid(col,row)

      snowmodel_struc(n)%sm(t)%lwdown = &
            Qli_grid(col,row)

      ! U-/V-wind components (units: m/s):
      snowmodel_struc(n)%sm(t)%uwind = &
             uwind_grid(col,row)
      snowmodel_struc(n)%sm(t)%vwind = &
             vwind_grid(col,row)

      ! Write out SnowModel metforcing file fields ("native" text file fields):
      if( snowmodel_struc(n)%write_sm_metfields == 1 ) then

        call LIS_diagnoseSurfaceOutputVar(n, t,LIS_MOC_TAIRFORC, &
             value=snowmodel_struc(n)%sm(t)%tair,&
             unit="K",vlevel=1,direction="-",&
             surface_type=LIS_rc%lsm_index)

        call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_SWDOWNFORC_SM,&
             value=snowmodel_struc(n)%sm(t)%swdown,&
             unit="W m-2",vlevel=1,direction="DN",&
             surface_type=LIS_rc%lsm_index)

        call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_LWDOWNFORC_SM,&
             value=snowmodel_struc(n)%sm(t)%lwdown,&
             unit="W m-2",vlevel=1,direction="DN",&
             surface_type=LIS_rc%lsm_index)

        call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_EWINDFORC_SM, &
             value=snowmodel_struc(n)%sm(t)%uwind, &
             unit="m s-1", vlevel=1, direction="E", &
             surface_type=LIS_rc%lsm_index)

        call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_NWINDFORC_SM, &
             value=snowmodel_struc(n)%sm(t)%vwind, &
             unit="m s-1", vlevel=1, direction="N", &
             surface_type=LIS_rc%lsm_index)
    
      endif

      ! -----------------------------------------

      ! SWE depth (units: m, SnowModel)
      snowmodel_struc(n)%sm(t)%swe_depth = &
            swe_depth(col,row)    

      call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_SWE_SM, &
           value=snowmodel_struc(n)%sm(t)%swe_depth, &
           unit="m", vlevel=1, direction="-", &
           surface_type=LIS_rc%lsm_index)

!      IF user wants to write out SWE in kg m-2 units (need conversion)
!      call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_SWE, &
!           value=snowmodel_struc(n)%sm(t)%swe_depth, &
!           unit="kg m-2", vlevel=1, direction="-", &
!           surface_type=LIS_rc%lsm_index)

      ! Snow depth (units: m) 
      snowmodel_struc(n)%sm(t)%snow_depth = &
            snow_depth(col,row)    
      call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_SNOWDEPTH_SM, &
           value=snowmodel_struc(n)%sm(t)%snow_depth, &
           unit="m", vlevel=1, direction="-", &
           surface_type = LIS_rc%lsm_index)

      ! Snow density (units: kg m-3):
      snowmodel_struc(n)%sm(t)%sden = &
            xro_snow(col,row)
      call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_SNOWDENSITY_SM, &
           value = snowmodel_struc(n)%sm(t)%sden, &
           unit="kg m-3", vlevel=1, direction="-", &
           surface_type = LIS_rc%lsm_index)

      ! SWE melt (units: m, for Snowmodel; units: kg m-2 s-1, for other models). 
      !  Snow melt (QSM):
      snowmodel_struc(n)%sm(t)%swemelt = &
            swemelt(col,row)
      call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_QSM_SM, &
           value = snowmodel_struc(n)%sm(t)%swemelt, &
           unit="m", vlevel=1, direction="S2L", &
!           unit="kg m-2 s-1", vlevel=1, direction="S2L", &
!          Note: "S2L" is "Solid to Liquid" 
!         (https://www.lmd.jussieu.fr/~polcher/ALMA/convention_output_3.html)
           surface_type = LIS_rc%lsm_index)

      ! Snow sublimation (units: m, for Snowmodel; kg m-2 s-1, for other models):
      snowmodel_struc(n)%sm(t)%sublim = &
            swesublim(col,row)
      call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_SUBSNOW_SM, &
           value = snowmodel_struc(n)%sm(t)%sublim, &
           unit="m", vlevel = 1, direction="-", &
!           unit="kg m-2 s-1", vlevel = 1, direction="-", &
           surface_type = LIS_rc%lsm_index )

      ! Runoff from base of snowpack 
      !  (units: m/dt for Snowmodel; otherwise, units: kg m-2 s-1):
      snowmodel_struc(n)%sm(t)%runoff = &
            runoff(col,row)
      call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_QS_SM, &
           value = snowmodel_struc(n)%sm(t)%runoff, &
           unit="m", vlevel=1, direction="OUT", &
!           unit="kg m-2 s-1", vlevel=1, direction="OUT", &
           surface_type = LIS_rc%lsm_index)

      ! Total (water-equivalent) precipitation (units: m/time_step, for SnowModel)
      snowmodel_struc(n)%sm(t)%totprecip = &
            prec_grid(col,row)
      call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_TOTALPRECIP_SM, &
           value=snowmodel_struc(n)%sm(t)%totprecip, &
           unit="m", vlevel=1, direction="DN", &
!           unit="kg m-2", vlevel=1, direction="DN", &
           surface_type=LIS_rc%lsm_index)

      snowmodel_struc(n)%sm(t)%totprecip = &
            prec_grid(col,row)
      call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_TOTALPRECIP_SM, &
           value=snowmodel_struc(n)%sm(t)%totprecip*1000., &
           unit="kg m-2", vlevel=1, direction="DN", &
           surface_type=LIS_rc%lsm_index)

      snowmodel_struc(n)%sm(t)%totprecip = &
            prec_grid(col,row)
      call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_TOTALPRECIP_SM, &
           value=snowmodel_struc(n)%sm(t)%totprecip*1000./dt, &
           unit="kg m-2 s-1", vlevel=1, direction="DN", &
           surface_type=LIS_rc%lsm_index)


      ! Liquid precipitation portion (units: m/time_step, for SnowModel)
      snowmodel_struc(n)%sm(t)%rainf = &
            rain(col,row)
      call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_RAINF_SM, &
           value=snowmodel_struc(n)%sm(t)%rainf, &
           unit="m", vlevel=1, direction="DN", &
!           unit="kg m-2", vlevel=1, direction="DN", &
           surface_type=LIS_rc%lsm_index)

      snowmodel_struc(n)%sm(t)%rainf = &
            rain(col,row)
      call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_RAINF_SM, &
           value=snowmodel_struc(n)%sm(t)%rainf*1000., &
           unit="kg m-2", vlevel=1, direction="DN", &
           surface_type=LIS_rc%lsm_index)

      snowmodel_struc(n)%sm(t)%rainf = &
            rain(col,row)
      call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_RAINF_SM, &
           value=snowmodel_struc(n)%sm(t)%rainf*1000./dt, &
           unit="kg m-2 s-1", vlevel=1, direction="DN", &
           surface_type=LIS_rc%lsm_index)


     ! Solid precipitation portion (units: m/time_step, for SnowModel)
     snowmodel_struc(n)%sm(t)%snowf = &
            sprec(col,row)
     call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_SNOWF_SM, &
           value=snowmodel_struc(n)%sm(t)%snowf, &
           unit="m", vlevel=1, direction="DN", &
!           unit="kg m-2", vlevel=1, direction="DN", &
           surface_type=LIS_rc%lsm_index)

      snowmodel_struc(n)%sm(t)%snowf = &
            sprec(col,row)
      call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_SNOWF_SM, &
           value=snowmodel_struc(n)%sm(t)%snowf*1000., &
           unit="kg m-2", vlevel=1, direction="DN", &
           surface_type=LIS_rc%lsm_index)

      snowmodel_struc(n)%sm(t)%snowf = &
            sprec(col,row)
      call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_SNOWF_SM, &
           value=snowmodel_struc(n)%sm(t)%snowf*1000./dt, &
           unit="kg m-2 s-1", vlevel=1, direction="DN", &
           surface_type=LIS_rc%lsm_index)


     ! Albedo (units: -)
     snowmodel_struc(n)%sm(t)%albedo = &
            albedo(col,row)
     call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_ALBEDO_SM, &
           value=snowmodel_struc(n)%sm(t)%albedo, &
           unit="-", vlevel=1, direction="-", &
           surface_type=LIS_rc%lsm_index)

     ! Parameter fields:
     ! Topographic - elevation (m):
     call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_ELEVATION_SM, &
!           value=snowmodel_struc(n)%sm(t)%smtopo, &
           value=topo(col,row), &
!           value=topo_land(col,row), &
!           value=curvature(col,row), &
!           value=slope_az(col,row), &
!           value=terrain_slope(col,row), &
           unit="m", vlevel=1, direction="-", &
           surface_type=LIS_rc%lsm_index)

     ! Vegetation type:
     call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_LANDCOVER_SM, &
!           value=snowmodel_struc(n)%sm(t)%smvege, &
           value=vegtype(col,row), &
           unit="-", vlevel=1, direction="-", &
           surface_type=LIS_rc%lsm_index)

! ..............................................................
! 
! Writing state values for writing out to restart files:

     snowmodel_struc(n)%sm(t)%snow_d = &
              snow_d(col,row)
     snowmodel_struc(n)%sm(t)%canopy_int = &
              canopy_int(col,row)
     snowmodel_struc(n)%sm(t)%soft_snow_d = &
              soft_snow_d(col,row)
     snowmodel_struc(n)%sm(t)%ro_snow_grid = &
              ro_snow_grid(col,row)
     snowmodel_struc(n)%sm(t)%ro_soft_snow_old = &
              ro_soft_snow_old(col,row)
     snowmodel_struc(n)%sm(t)%snow_d_init = &
              snow_d_init(col,row)
     snowmodel_struc(n)%sm(t)%swe_depth_old = &
              swe_depth_old(col,row)
     snowmodel_struc(n)%sm(t)%canopy_int_old = &
              canopy_int_old(col,row)
     snowmodel_struc(n)%sm(t)%topo = &
              topo(col,row)
     snowmodel_struc(n)%sm(t)%sum_sprec = &
              sum_sprec(col,row)

!
! FIELDS THAT LIS CAN SUPPORT AND MAYBE AVAILABLE BY SNOWMODEL
!  BUT ARE NOT WRITTEN OUT ...
      ! Latent hear flux 
!   snowmodel_struc(n)%sm%qe = 40.0
!      call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_QLE, &
!           value=snowmodel_struc(n)%sm(t)%qe, &
!           vlevel=1,unit="W m-2",direction="UP",surface_type=LIS_rc%lsm_index)

      ! Snowliq (unit=mm ). ***  snow-layer liquid water 
!      do i=1, snowmodel_struc(n)%nsnow
!         call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_SNOWLIQ, &
!              value = NOAHMP36_struc(n)%noahmp36(t)%snowliq(i), &
!              unit="mm", vlevel=i, direction="-", &
!              surface_type = LIS_rc%lsm_index )
!      enddo

     ! Snowliq (unit=mm ). ***  snow-layer ice
!      do i=1, snowmodel_struc(n)%nsnow
!         call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_SNOWICE, &
!              value = NOAHMP36_struc(n)%noahmp36(t)%snowce(i), &
!              unit="mm", vlevel=i, direction="-", &
!              surface_type = LIS_rc%lsm_index )
!      enddo

     ! Snow_temp (unit=K). ***  snow layer temperature
!      do i=1, NOAHMP36_struc(n)%nsnow
! Test code to reset snow temperature to undefined
! when there is no corresponding snow layer - Mocko
!     if ((i + abs(NOAHMP36_struc(n)%noahmp36(t)%isnow))     &
!        .le.NOAHMP36_struc(n)%nsnow) then
!         snow_temp(i) = LIS_rc%udef
!     endif
!     call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_SNOWTPROF, &
!          value = snow_temp(i), &
!          unit="K", vlevel=i, direction="-", &
!          surface_type = LIS_rc%lsm_index)
!      end do

      ! Snow density for each layer
!     !  ro_layer or sdenz .. 
!      do i=1, NOAHMP36_struc(n)%nsnow
!         call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_LAYERSNOWDENSITY, 
!         value = bdsno,&
!         vlevel=i, unit="kg m-3", direction="-", &
!         surface_type = LIS_rc%lsm_index)
!      enddo
!
! .........................................................................

       ! Reset forcing variables to 0 and counter:
       snowmodel_struc(n)%forc_count = 0
       snowmodel_struc(n)%sm(t)%tair = 0
       snowmodel_struc(n)%sm(t)%qair = 0
       snowmodel_struc(n)%sm(t)%swdown = 0
       snowmodel_struc(n)%sm(t)%lwdown = 0
       snowmodel_struc(n)%sm(t)%uwind = 0
       snowmodel_struc(n)%sm(t)%vwind = 0
       snowmodel_struc(n)%sm(t)%psurf = 0
       snowmodel_struc(n)%sm(t)%rainf = 0
       snowmodel_struc(n)%sm(t)%snowf = 0

   enddo   ! End 1-d tile space loop

   ! For multi-year simulations, sometimes it is desirable to zero
   !   out the snow cover arrays on a certain summer date, to prevent
   !   glaciers from forming.

   if (imonth.eq.iclear_mn .and. iday.eq.iclear_dy .and.&
        xhour.eq.xclear_hr) then

      if (seaice_run.eq.4.0) then
         CALL ZERO_SNOW_SEAICE_4(LIS_rc%lnc(n),LIS_rc%lnr(n),snow_depth,ro_snow_grid,&
            ro_snow,swe_depth,swe_depth_old,canopy_int_old,KK,&
            sum_swemelt,tslsnowfall,snod_layer,swed_layer,&
            ro_layer,T_old,sum_sprec,multilayer_snowpack,&
            tsls_threshold,iyear,diam_layer,output_path_wo_assim)
      else
         CALL ZERO_SNOW(LIS_rc%lnc(n),LIS_rc%lnr(n),snow_depth,ro_snow_grid,ro_snow,&
            swe_depth,swe_depth_old,canopy_int_old,KK,sum_swemelt,&
            tslsnowfall,snod_layer,swed_layer,ro_layer,T_old,&
            sum_sprec,multilayer_snowpack,tsls_threshold,&
            sum_trans)
      endif

   endif

   if( iter == max_iter ) then
     write(LIS_logunit,*) "[INFO] Reached max iteration count: ",iter
   endif

end subroutine snowmodel_main
