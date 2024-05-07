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
! !ROUTINE: snowmodel_setup
! \label{snowmodel_setup}
!
! !REVISION HISTORY:
!  14 Apr 2020: Kristi Arsenault; Add G. Liston's SnowModel 
! 
! !INTERFACE:
subroutine snowmodel_setup()
! !USES:
   use LIS_logMod,    only : LIS_verify, LIS_logunit, LIS_endrun
   use LIS_coreMod,   only : LIS_rc, LIS_surface,&
         LIS_localPet, LIS_ews_ind, LIS_ewe_ind,&
         LIS_nss_ind, LIS_nse_ind, LIS_ews_halo_ind, LIS_ewe_halo_ind,&
         LIS_nss_halo_ind, LIS_nse_halo_ind
   use LIS_fileIOMod, only : LIS_read_param
   use snowmodel_lsmMod
   use snowmodel_module
   use snowmodel_inc
   use snowmodel_vars
!
! !DESCRIPTION: 
! 
!  This routine is the entry point to set up the parameters
!  required for Liston's SnowModel. These include vegetation
!  topography, and initialization of state variables in SnowModel.
!  
! The routines invoked are: 
! \begin{description}
! \item[snowmodel\_setvegparms](\ref{snowmodel_setvegparms}) \newline
!   initializes the vegetation-related parameters in SnowModel
! \end{description}
!EOP

  implicit none
  integer           :: n
  integer           :: mtype
  integer           :: t
  integer           :: col, row
  integer           :: ews, ewe, nss, nse 
  double precision  :: xmn_part  ! center x of local LL starting point
  double precision  :: ymn_part  ! center y of local LL starting point
  integer, allocatable :: global_kstn(:,:,:)
  real, allocatable :: placeholder(:,:)

! _______________________________________________________________

  mtype = LIS_rc%lsm_index

  write(LIS_logunit,*)"[INFO] Reading in SnowModel vegetation and topography (snowmodel_setup)"

  ! Read in spatial SnowModel maps from LDT:
  do n=1, LIS_rc%nnest

     ! Allocate memory for place holder for #n nest
     allocate(placeholder(LIS_rc%lnc(n), LIS_rc%lnr(n)))

     ews = LIS_ews_halo_ind(n, LIS_localPet+1) 
     ewe = LIS_ewe_halo_ind(n, LIS_localPet+1)
     nss = LIS_nss_halo_ind(n, LIS_localPet+1)
     nse = LIS_nse_halo_ind(n, LIS_localPet+1)

     ! Initialize topography and landcover parameters:
     snowmodel_struc(n)%sm(:)%smtopo = 0.
     snowmodel_struc(n)%sm(:)%smvege = 0.

     ! Assign LDT-based landcover and topographic maps :
     if( snowmodel_struc(n)%sm_params_opt == "LDT" ) then
        write(LIS_logunit,*) "[INFO] Reading in SnowModel LSM parameters from LDT"
        ascii_topoveg = 2.0    ! No file read in by SnowModel; use LDT input
     else
        write(LIS_logunit,*) "[INFO] Reading in SnowModel LSM parameters from snowmodel.par file "
     endif

     if( ascii_topoveg == 2.0 ) then
       write(LIS_logunit,*) "SnowModel: Reading parameter SMTOPO from: ",&
          trim(LIS_rc%paramfile(n))
       call LIS_read_param(n, "SMTOPO", placeholder)
       do t = 1, LIS_rc%npatch(n, mtype)
          col = LIS_surface(n, mtype)%tile(t)%col
          row = LIS_surface(n, mtype)%tile(t)%row
          snowmodel_struc(n)%sm(t)%smtopo = placeholder(col, row)
         ! Assign LDT-version of SnowModel topo map to topo_land ...
          topo_land(col,row) = placeholder(col, row)
       enddo

       write(LIS_logunit,*) "SnowModel: Reading parameter SMVEG from: ",&
          trim(LIS_rc%paramfile(n))
       call LIS_read_param(n, "SMVEG", placeholder)
       do t = 1, LIS_rc%npatch(n, mtype)
          col = LIS_surface(n, mtype)%tile(t)%col
          row = LIS_surface(n, mtype)%tile(t)%row
          snowmodel_struc(n)%sm(t)%smvege = placeholder(col, row)
         ! Assign LDT-version of SnowModel veg map to vegtype ...
          vegtype(col,row) = placeholder(col, row)
       enddo

     endif

  enddo

!  call snowmodel_setvegparms(LIS_rc%lsm_index)

  ! ------------------------------

  ! SnowModel Main code calls and options for reading in and 
  !  preprocessing necessary inputs to the submodel components:
  ! All below from the original snowmodel_main.f ...

  do n=1,LIS_rc%nnest

     write(LIS_logunit,*) "[INFO] SnowModel 'main' calls for preprocessing inputs"

     ! Using local parallel subdomain starting index values:
     ! LIS_ews_halo_ind(n,LIS_localPet+1) -- defined in LIS_coreMod.F90 ...
     xmn_part = xmn + deltax * ( real(LIS_ews_halo_ind(n,LIS_localPet+1)) - 1.0 )
     ymn_part = ymn + deltay * ( real(LIS_nss_halo_ind(n,LIS_localPet+1)) - 1.0 )

     ! This loop runs the correction/data assimilation adjustment
     !  iterations.
     if (ihrestart_flag.ge.0) then
       if (i_dataassim_loop.lt.0.0) then
         i_corr_start = 2
       else
         i_corr_start = 1
       endif
     else
       i_corr_start = 1
     endif

     do icorr_factor_loop=i_corr_start,irun_data_assim+1

       ! Perform a variety of preprocessing and model setup steps, like
       !  read in topography and vegetation arrays, open input and output
       !  files, etc.

       if( snowmodel_struc(n)%call_sm_preproc == 1 ) then

         write(LIS_logunit,*) "[INFO] Calling original 'preprocess_code' routine "

         CALL PREPROCESS_CODE(topoveg_fname,const_veg_flag,&
          vegtype,veg_z0,vegsnowdepth,fetch,xmu,C_z,h_const,&   
          wind_min,Up_const,dz_susp,ztop_susp,fall_vel,Ur_const,&
          ro_water,ro_air,gravity,vonKarman,pi,twopio360,snow_z0,&
          LIS_rc%lnc(n),LIS_rc%lnr(n),sum_sprec,sum_qsubl,sum_trans,sum_unload,topo,&  ! KRA: nx,ny
          topo_land,snow_d,topoflag,snow_d_init,snow_d_init_const,&  
          soft_snow_d,met_input_fname,igrads_metfile,deltax,deltay,&
          snowtran_output_fname,micromet_output_fname,&
          enbal_output_fname,snowpack_output_fname,print_micromet,&
          print_enbal,print_snowpack,print_snowtran,run_micromet,&
          run_enbal,run_snowpack,run_snowtran,ro_snow_grid,swe_depth,&
          sum_runoff,sum_prec,ro_snow,twolayer_flag,sum_Qcs,&
          canopy_int,ascii_topoveg,topo_ascii_fname,icorr_factor_loop,&
          veg_ascii_fname,undef,isingle_stn_flag,max_iter,&
          i_tair_flag,i_rh_flag,i_wind_flag,i_prec_flag,sum_glacmelt,&
          snow_depth,sum_d_canopy_int,corr_factor,icorr_factor_index,&
          sum_sfcsublim,barnes_lg_domain,n_stns_used,k_stn,xmn_part,ymn_part,& !KRA: xmn, ymn
          ro_soft_snow_old,sum_swemelt,xlat,lat_solar_flag,xlat_grid,&
          xlon_grid,UTC_flag,dt,swe_depth_old,canopy_int_old,&
          vegsnowd_xy,iveg_ht_flag,ihrestart_flag,i_dataassim_loop,&
          multilayer_snowpack,max_layers,multilayer_output_fname,&
          print_multilayer,KK,tslsnowfall,tsls_threshold,&
          irun_data_assim,izero_snow_date,iclear_mn,iclear_dy,&
          xclear_hr,snod_layer,swed_layer,ro_layer,T_old,gamma,&
          icond_flag,curve_lg_scale_flag,curve_wt_lg,check_met_data,&
          seaice_run,snowmodel_line_flag,xg_line,yg_line,print_user,&
          cf_precip_flag,cf_precip,print_inc,xhour_init,Tabler_1_flag,&
          Tabler_2_flag,iyear_init,imonth_init,iday_init,print_var,&
          output_path_wo_assim,output_path_wi_assim,nrecs_max,&
          tabler_sfc_path_name,print_outvars,diam_layer)


        ! Generate kstn arrays for MicroMet inputs (on SnowModel side):
        if( snowmodel_struc(n)%sm_micromet_opt == "SnowModel" ) then

          allocate( global_kstn(LIS_rc%gnc(n),LIS_rc%gnr(n),9) )
          global_kstn = 0

          ! If the large-domain barnes oi scheme is used, generate the
          !   nearest-station indexing array.
          if (barnes_lg_domain.eq.1.0) then
            if (n_stns_used.gt.9 .or. n_stns_used.lt.1) then
              print *,'invalid n_stns_used value'
              stop
            endif
! Original call from preprocess.f (which is commented now in that routine):
!          call get_nearest_stns_1(nx,ny,xmn,ymn,deltax,deltay, &
!                n_stns_used,k_stn,snowmodel_line_flag,xg_line,yg_line)
! Updated here to represent entire SnowModel run domain:
            call get_nearest_stns_1( LIS_rc%gnc(n),LIS_rc%gnr(n),&
                  xmn, ymn, deltax, deltay,&
                  n_stns_used, global_kstn, &
                  snowmodel_line_flag, xg_line, yg_line)

            k_stn(:,:,:) = global_kstn(ews:ewe,nss:nse,:)
          endif
          deallocate(global_kstn)

        elseif( snowmodel_struc(n)%sm_micromet_opt == "LIS" ) then
            call get_nearest_stns_1( LIS_rc%lnc(n),LIS_rc%lnr(n),&
                  xmn_part, ymn_part, deltax, deltay,&
                  n_stns_used, k_stn, &
                  snowmodel_line_flag, xg_line, yg_line)

        endif  ! MM option: Snowmodel or LIS
       endif
     end do

   end do   ! Nest loop

 
end subroutine snowmodel_setup
 
