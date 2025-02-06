! readparam_code.f

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine READPARAM_CODE(dt,deltax,deltay,Utau_t_flag,&
     &  subgrid_flag,twolayer_flag,snowmodel_dot_par_fname,&
     &  bc_flag,curve_len_scale,slopewt,curvewt,ht_windobs,&
     &  ht_rhobs,ro_snow,snow_d_init_const,const_veg_flag,&
     &  vegsnowdepth,nx,ny,max_iter,met_input_fname,xmn,ymn,&
     &  iyear_init,imonth_init,iday_init,xhour_init,undef,ifill,&
     &  iobsint,dn,xlat,i_tair_flag,i_rh_flag,i_wind_flag,&
     &  i_solar_flag,i_prec_flag,isingle_stn_flag,igrads_metfile,&
     &  windspd_min,icond_flag,run_micromet,run_enbal,run_snowpack,&
     &  run_snowtran,topoflag,topoveg_fname,snowtran_output_fname,&
     &  micromet_output_fname,enbal_output_fname,Utau_t_const,&
     &  snowpack_output_fname,print_micromet,print_enbal,&
     &  print_snowpack,print_snowtran,i_longwave_flag,print_user,&
     &  ascii_topoveg,topo_ascii_fname,veg_ascii_fname,&
     &  irun_data_assim,lapse_rate_user_flag,&
     &  iprecip_lapse_rate_user_flag,use_shortwave_obs,&
     &  use_longwave_obs,use_sfc_pressure_obs,calc_subcanopy_met,&
     &  sfc_sublim_flag,gap_frac,cloud_frac_factor,&
     &  albedo_snow_forest,albedo_snow_clearing,albedo_glacier,&
     &  barnes_lg_domain,n_stns_used,tabler_dir,slope_adjust,&
     &  lat_solar_flag,UTC_flag,iveg_ht_flag,ihrestart_flag,&
     &  ihrestart_inc,i_dataassim_loop,tsls_threshold,dz_snow_min,&
     &  print_multilayer,multilayer_snowpack,max_layers,&
     &  multilayer_output_fname,izero_snow_date,curve_lg_scale_flag,&
     &  check_met_data,seaice_run,snowmodel_line_flag,wind_lapse_rate,&
     &  iprecip_scheme,cf_precip_flag,snowfall_frac,print_inc,&
     &  output_path_wo_assim,output_path_wi_assim,Tabler_1_flag,&
     &  Tabler_2_flag,tabler_sfc_path_name,print_var)

! These programs read and process the input parameter data.
!
! The following must be true:
!
!   All comment lines start with a ! in the first position.
!
!   Blank lines are permitted.
!
!   All parameter statements start with the parameter name, followed
!   by a space, followed by an =, followed by a space, followed by
!   the actual value, with nothing after that.  These statements can
!   have leading blanks, and must fit within 120 columns, and the
!   parameter value cannot be longer than 80 characters (so, for
!   example, your paths/file_names to the right of the "=" sign
!   cannot be longer than 80 characters).
!
!   It is set up to deal with integer, real, and string (character)
!   input values.

      use snowmodel_vars, only: allocate_arrays
      use snowmodel_inc
      implicit none

! Put parameter names here:
      real dt,deltax,deltay,Utau_t_flag,topoflag,Utau_t_const,&
     &  subgrid_flag,twolayer_flag,&
     &  bc_flag,curve_len_scale,slopewt,curvewt,ht_windobs,&
     &  ht_rhobs,ro_snow,snow_d_init_const,const_veg_flag,&
     &  windspd_min,ascii_topoveg,gap_frac,cloud_frac_factor,&
     &  albedo_snow_forest,albedo_snow_clearing,albedo_glacier,&
     &  barnes_lg_domain,tabler_dir,slope_adjust,UTC_flag,&
     &  tsls_threshold,dz_snow_min,print_multilayer,&
     &  curve_lg_scale_flag,check_met_data,seaice_run,&
     &  snowmodel_line_flag,wind_lapse_rate,cf_precip_flag,&
     &  snowfall_frac,Tabler_1_flag,Tabler_2_flag

      character*1 print_var(n_print_vars)

      character*1 print_var_01,print_var_02,print_var_03,&
     &  print_var_04,print_var_05,print_var_06,print_var_07,&
     &  print_var_08,print_var_09,print_var_10,print_var_11,&
     &  print_var_12,print_var_13,print_var_14,print_var_15,&
     &  print_var_16,print_var_17,print_var_18,print_var_19,&
     &  print_var_20,print_var_21,print_var_22,print_var_23,&
     &  print_var_24,print_var_25,print_var_26,print_var_27,&
     &  print_var_28,print_var_29,print_var_30

      real vegsnowdepth(nvegtypes)
      real run_micromet,run_enbal,run_snowpack,run_snowtran,&
     &  print_micromet,print_enbal,print_snowpack,print_snowtran,&
     &  print_user,use_shortwave_obs,use_longwave_obs,print_inc,&
     &  use_sfc_pressure_obs,calc_subcanopy_met,sfc_sublim_flag

      integer nx,ny,max_iter,icond_flag,irun_data_assim

      character*80 topoveg_fname,met_input_fname,topo_ascii_fname,&
     &  veg_ascii_fname

      character*80 snowtran_output_fname,micromet_output_fname,&
     &  enbal_output_fname,snowpack_output_fname,multilayer_output_fname

      character*80 output_path_wo_assim,output_path_wi_assim,&
     &  tabler_sfc_path_name

      integer iyear_init,imonth_init,iday_init,&
     &  i_tair_flag,i_rh_flag,i_wind_flag,i_solar_flag,i_prec_flag,&
     &  i_longwave_flag,isingle_stn_flag,igrads_metfile,&
     &  lapse_rate_user_flag,iprecip_lapse_rate_user_flag,&
     &  n_stns_used,iveg_ht_flag,lat_solar_flag,ihrestart_inc,&
     &  ihrestart_flag,i_dataassim_loop,multilayer_snowpack,max_layers,&
     &  izero_snow_date,iprecip_scheme, metforce_opt

      double precision xmn,ymn
      real xhour_init,dn
      real undef               ! undefined value
      integer ifill    ! flag (=1) forces a value in every cell
      integer iobsint  ! flag (=1) use dn value from .par file
      real xlat      ! approx. latitude of domain center, decimal deg

      character*100 snowmodel_dot_par_fname

! Working parameters.
      character*120 input_string
      character*80 c_param
      character*80 c_value
      integer k,max_par_lines,i_param_chars,i_value_chars,&
     &  icomment_flag,ipar_flag

!      integer, parameter :: npars=139
      integer, parameter :: npars=140
      integer ipar_count
      character*40 cpar_name(npars)

      max_par_lines = 1500

!     open (40,file='snowmodel.par')
      open (40,file=snowmodel_dot_par_fname)

! Initialize the input-parameter counting array.  This is used to
!   make sure all of the input parameters are read in correctly.
!   Also initialize a array that will collect the variable names
!   that have been read in, so it can be compared against the .par
!   file to help figure out which variables have not been defined
!   before the model simulation starts.
      ipar_count = 0
      do k=1,npars
        cpar_name(k) = '  BLANK: SOMETHING IS NOT DEFINED'
      enddo

      do k=1,max_par_lines
        read (40,'(a120)',end=99) input_string

        call get_param_data(input_string,c_param,c_value,&
     &    i_param_chars,i_value_chars,icomment_flag)

! Process the string if it is not a comment.
        if (icomment_flag.eq.0) then

! GENERAL MODEL SETUP.
          if (c_param(1:i_param_chars).eq.'nx') then
            ipar_count = ipar_count + 1
            cpar_name(ipar_count) = c_param(1:i_param_chars)
! KRA
!            call char2int(nx,i_value_chars,c_value,&
!     &        c_param(1:i_param_chars))
            if (nx.le.0) then
              print *,'nx cannot be less than or equal to 0'
              stop
            endif
          endif

          if (c_param(1:i_param_chars).eq.'ny') then
            ipar_count = ipar_count + 1
            cpar_name(ipar_count) = c_param(1:i_param_chars)
! KRA
!            call char2int(ny,i_value_chars,c_value,&
!     &        c_param(1:i_param_chars))
            if (ny.le.0) then
              print *,'ny cannot be less than or equal to 0'
              stop
            endif
          endif

          if (c_param(1:i_param_chars).eq.'deltax') then
            ipar_count = ipar_count + 1
            cpar_name(ipar_count) = c_param(1:i_param_chars)
            call char2real(deltax,i_value_chars,c_value,&
     &        c_param(1:i_param_chars))
            if (deltax.le.0.0) then
              print *,'deltax cannot be less than or equal to 0.0'
              stop
            endif
          endif

          !Subroutine effective only on first invocation.
          !Check performed inside the routine.
          if(nx /= 0 .and. ny /= 0) call allocate_arrays(nx,ny)

          if (c_param(1:i_param_chars).eq.'deltay') then
            ipar_count = ipar_count + 1
            cpar_name(ipar_count) = c_param(1:i_param_chars)
            call char2real(deltay,i_value_chars,c_value,&
     &        c_param(1:i_param_chars))
            if (deltay.le.0.0) then
              print *,'deltay cannot be less than or equal to 0.0'
              stop
            endif
          endif

          if (c_param(1:i_param_chars).eq.'xmn') then
            ipar_count = ipar_count + 1
            cpar_name(ipar_count) = c_param(1:i_param_chars)
            call char2double(xmn,i_value_chars,c_value,&
     &        c_param(1:i_param_chars))
          endif

          if (c_param(1:i_param_chars).eq.'ymn') then
            ipar_count = ipar_count + 1
            cpar_name(ipar_count) = c_param(1:i_param_chars)
            call char2double(ymn,i_value_chars,c_value,&
     &        c_param(1:i_param_chars))
          endif

          if (c_param(1:i_param_chars).eq.'dt') then
            ipar_count = ipar_count + 1
            cpar_name(ipar_count) = c_param(1:i_param_chars)
            call char2real(dt,i_value_chars,c_value,&
     &        c_param(1:i_param_chars))
            if (dt.le.0.0 .or. dt.gt.86400.0) then
              print *,'dt out of range'
              stop
            endif
          endif

          if (c_param(1:i_param_chars).eq.'iyear_init') then
            ipar_count = ipar_count + 1
            cpar_name(ipar_count) = c_param(1:i_param_chars)
            call char2int(iyear_init,i_value_chars,c_value,&
     &        c_param(1:i_param_chars))
            if (iyear_init.lt.1900 .or. iyear_init.gt.2100) then
              print *,'iyear_init out of range'
              stop
            endif
          endif

          if (c_param(1:i_param_chars).eq.'imonth_init') then
            ipar_count = ipar_count + 1
            cpar_name(ipar_count) = c_param(1:i_param_chars)
            call char2int(imonth_init,i_value_chars,c_value,&
     &        c_param(1:i_param_chars))
            if (imonth_init.lt.1 .or. imonth_init.gt.12) then
              print *,'imonth_init out of range'
              stop
            endif
          endif

          if (c_param(1:i_param_chars).eq.'iday_init') then
            ipar_count = ipar_count + 1
            cpar_name(ipar_count) = c_param(1:i_param_chars)
            call char2int(iday_init,i_value_chars,c_value,&
     &        c_param(1:i_param_chars))
            if (iday_init.lt.1 .or. iday_init.gt.31) then
              print *,'iday_init out of range'
              stop
            endif
          endif

          if (c_param(1:i_param_chars).eq.'xhour_init') then
            ipar_count = ipar_count + 1
            cpar_name(ipar_count) = c_param(1:i_param_chars)
            call char2real(xhour_init,i_value_chars,c_value,&
     &        c_param(1:i_param_chars))
            if (xhour_init.lt.0.0 .or. xhour_init.ge.24.0) then
              print *,'xhour_init out of range'
              stop
            endif
          endif

          if (c_param(1:i_param_chars).eq.'max_iter') then
            ipar_count = ipar_count + 1
            cpar_name(ipar_count) = c_param(1:i_param_chars)
            call char2int(max_iter,i_value_chars,c_value,&
     &        c_param(1:i_param_chars))
            if (max_iter.le.0 .or. max_iter.gt.300000) then
              print *,'max_iter out of range'
              stop
            endif
          endif

          if (c_param(1:i_param_chars).eq.'isingle_stn_flag') then
            ipar_count = ipar_count + 1
            cpar_name(ipar_count) = c_param(1:i_param_chars)
            call char2int(isingle_stn_flag,i_value_chars,c_value,&
     &        c_param(1:i_param_chars))
            if (isingle_stn_flag.ne.0 .and. isingle_stn_flag.ne.1) then
              print *,'isingle_stn_flag not 0 or 1'
              stop
            endif
          endif

          if (c_param(1:i_param_chars).eq.'metforce_opt') then
            ipar_count = ipar_count + 1
            cpar_name(ipar_count) = c_param(1:i_param_chars)
            call char2int(metforce_opt,i_value_chars,c_value,&
     &        c_param(1:i_param_chars))
            if (metforce_opt.ne.0 .and. metforce_opt.ne.1) then
              print *,'metforce_opt not 0 or 1'
              stop
            endif
          endif


          if (c_param(1:i_param_chars).eq.'igrads_metfile') then
            ipar_count = ipar_count + 1
            cpar_name(ipar_count) = c_param(1:i_param_chars)
            call char2int(igrads_metfile,i_value_chars,c_value,&
     &        c_param(1:i_param_chars))
            if (igrads_metfile.ne.0 .and. igrads_metfile.ne.1) then
              print *,'igrads_metfile not 0 or 1'
              stop
            endif
          endif

          if (c_param(1:i_param_chars).eq.'met_input_fname') then
            ipar_count = ipar_count + 1
            cpar_name(ipar_count) = c_param(1:i_param_chars)
            call char2char(met_input_fname,c_value,i_value_chars,&
     &        c_param(1:i_param_chars))
          endif

          if (c_param(1:i_param_chars).eq.'undef') then
            ipar_count = ipar_count + 1
            cpar_name(ipar_count) = c_param(1:i_param_chars)
            call char2real(undef,i_value_chars,c_value,&
     &        c_param(1:i_param_chars))
          endif

          if (c_param(1:i_param_chars).eq.'ascii_topoveg') then
            ipar_count = ipar_count + 1
            cpar_name(ipar_count) = c_param(1:i_param_chars)
            call char2real(ascii_topoveg,i_value_chars,c_value,&
     &        c_param(1:i_param_chars))
! KRA
!           Update here to allow for fields to be read in from LDT (KRA; 2020-07-16)
            if (ascii_topoveg.lt.0.0 .and. ascii_topoveg.gt.2.0) then
               print *,'ascii_topoveg not 0.0, 1.0, or 2.0'
!            if (ascii_topoveg.ne.0.0 .and. ascii_topoveg.ne.1.0) then
!              print *,'ascii_topoveg not 0.0 or 1.0'
! KRA
              stop
            endif
          endif

          if (c_param(1:i_param_chars).eq.'topoveg_fname') then
            ipar_count = ipar_count + 1
            cpar_name(ipar_count) = c_param(1:i_param_chars)
            call char2char(topoveg_fname,c_value,i_value_chars,&
     &        c_param(1:i_param_chars))
          endif

          if (c_param(1:i_param_chars).eq.'topo_ascii_fname') then
            ipar_count = ipar_count + 1
            cpar_name(ipar_count) = c_param(1:i_param_chars)
            call char2char(topo_ascii_fname,c_value,i_value_chars,&
     &        c_param(1:i_param_chars))
          endif

          if (c_param(1:i_param_chars).eq.'veg_ascii_fname') then
            ipar_count = ipar_count + 1
            cpar_name(ipar_count) = c_param(1:i_param_chars)
            call char2char(veg_ascii_fname,c_value,i_value_chars,&
     &        c_param(1:i_param_chars))
          endif

!         if (c_param(1:i_param_chars).eq.'veg_shd_24') then
!           call char2real(vegsnowdepth(24),i_value_chars,c_value,&
!    &        c_param(1:i_param_chars))
!           if (vegsnowdepth(24).lt.0.0 .or. vegsnowdepth(24).gt.20.0)
!    &        then
!             print *,'veg_shd_24 out of range'
!             stop
!           endif
!         endif

          if (c_param(1:i_param_chars).eq.'veg_shd_25') then
            ipar_count = ipar_count + 1
            cpar_name(ipar_count) = c_param(1:i_param_chars)
            call char2real(vegsnowdepth(25),i_value_chars,c_value,&
     &        c_param(1:i_param_chars))
            if (vegsnowdepth(25).lt.0.0 .or. vegsnowdepth(25).gt.20.0)&
     &        then
              print *,'veg_shd_25 out of range'
              stop
            endif
          endif

          if (c_param(1:i_param_chars).eq.'veg_shd_26') then
            ipar_count = ipar_count + 1
            cpar_name(ipar_count) = c_param(1:i_param_chars)
            call char2real(vegsnowdepth(26),i_value_chars,c_value,&
     &        c_param(1:i_param_chars))
            if (vegsnowdepth(26).lt.0.0 .or. vegsnowdepth(26).gt.20.0)&
     &        then
              print *,'veg_shd_26 out of range'
              stop
            endif
          endif

          if (c_param(1:i_param_chars).eq.'veg_shd_27') then
            ipar_count = ipar_count + 1
            cpar_name(ipar_count) = c_param(1:i_param_chars)
            call char2real(vegsnowdepth(27),i_value_chars,c_value,&
     &        c_param(1:i_param_chars))
            if (vegsnowdepth(27).lt.0.0 .or. vegsnowdepth(27).gt.20.0)&
     &        then
              print *,'veg_shd_27 out of range'
              stop
            endif
          endif

          if (c_param(1:i_param_chars).eq.'veg_shd_28') then
            ipar_count = ipar_count + 1
            cpar_name(ipar_count) = c_param(1:i_param_chars)
            call char2real(vegsnowdepth(28),i_value_chars,c_value,&
     &        c_param(1:i_param_chars))
            if (vegsnowdepth(28).lt.0.0 .or. vegsnowdepth(28).gt.20.0)&
     &        then
              print *,'veg_shd_28 out of range'
              stop
            endif
          endif

          if (c_param(1:i_param_chars).eq.'veg_shd_29') then
            ipar_count = ipar_count + 1
            cpar_name(ipar_count) = c_param(1:i_param_chars)
            call char2real(vegsnowdepth(29),i_value_chars,c_value,&
     &        c_param(1:i_param_chars))
            if (vegsnowdepth(29).lt.0.0 .or. vegsnowdepth(29).gt.20.0)&
     &        then
              print *,'veg_shd_29 out of range'
              stop
            endif
          endif

          if (c_param(1:i_param_chars).eq.'veg_shd_30') then
            ipar_count = ipar_count + 1
            cpar_name(ipar_count) = c_param(1:i_param_chars)
            call char2real(vegsnowdepth(30),i_value_chars,c_value,&
     &        c_param(1:i_param_chars))
            if (vegsnowdepth(30).lt.0.0 .or. vegsnowdepth(30).gt.20.0)&
     &        then
              print *,'veg_shd_30 out of range'
              stop
            endif
          endif

          if (c_param(1:i_param_chars).eq.'const_veg_flag') then
            ipar_count = ipar_count + 1
            cpar_name(ipar_count) = c_param(1:i_param_chars)
            call char2real(const_veg_flag,i_value_chars,c_value,&
     &        c_param(1:i_param_chars))
            if (const_veg_flag.lt.0.0 .or. const_veg_flag.gt.30.0) then
              print *,'const_veg_flag out of range'
              stop
            endif
          endif

          if (c_param(1:i_param_chars).eq.'iveg_ht_flag') then
            ipar_count = ipar_count + 1
            cpar_name(ipar_count) = c_param(1:i_param_chars)
            call char2int(iveg_ht_flag,i_value_chars,c_value,&
     &        c_param(1:i_param_chars))
            if (iveg_ht_flag.lt.-1 .or. iveg_ht_flag.gt.1) then
              print *,'iveg_ht_flag out of range'
              stop
            endif
          endif

          if (c_param(1:i_param_chars).eq.'xlat') then
            ipar_count = ipar_count + 1
            cpar_name(ipar_count) = c_param(1:i_param_chars)
            call char2real(xlat,i_value_chars,c_value,&
     &        c_param(1:i_param_chars))
            if (xlat.lt.-90.0 .or. xlat.gt.90.0) then
              print *,'xlat out of range'
              stop
            endif
          endif

          if (c_param(1:i_param_chars).eq.'lat_solar_flag') then
            ipar_count = ipar_count + 1
            cpar_name(ipar_count) = c_param(1:i_param_chars)
            call char2int(lat_solar_flag,i_value_chars,c_value,&
     &        c_param(1:i_param_chars))
            if (lat_solar_flag.lt.-1 .or. lat_solar_flag.gt.1) then
              print *,'lat_solar_flag out of range'
              stop
            endif
          endif

          if (c_param(1:i_param_chars).eq.'UTC_flag') then
            ipar_count = ipar_count + 1
            cpar_name(ipar_count) = c_param(1:i_param_chars)
            call char2real(UTC_flag,i_value_chars,c_value,&
     &        c_param(1:i_param_chars))
            if (UTC_flag.lt.-1 .or. UTC_flag.gt.1) then
              print *,'UTC_flag out of range'
              stop
            endif
          endif

          if (c_param(1:i_param_chars).eq.'run_micromet') then
            ipar_count = ipar_count + 1
            cpar_name(ipar_count) = c_param(1:i_param_chars)
            call char2real(run_micromet,i_value_chars,c_value,&
     &        c_param(1:i_param_chars))
            if (run_micromet.ne.0.0 .and. run_micromet.ne.1.0) then
              print *,'run_micromet not 0.0 or 1.0'
              stop
            endif
          endif

          if (c_param(1:i_param_chars).eq.'run_enbal') then
            ipar_count = ipar_count + 1
            cpar_name(ipar_count) = c_param(1:i_param_chars)
            call char2real(run_enbal,i_value_chars,c_value,&
     &        c_param(1:i_param_chars))
            if (run_enbal.ne.0.0 .and. run_enbal.ne.1.0) then
              print *,'run_enbal not 0.0 or 1.0'
              stop
            endif
          endif

          if (c_param(1:i_param_chars).eq.'run_snowpack') then
            ipar_count = ipar_count + 1
            cpar_name(ipar_count) = c_param(1:i_param_chars)
            call char2real(run_snowpack,i_value_chars,c_value,&
     &        c_param(1:i_param_chars))
            if (run_snowpack.ne.0.0 .and. run_snowpack.ne.1.0) then
              print *,'run_snowpack not 0.0 or 1.0'
              stop
            endif
          endif

          if (c_param(1:i_param_chars).eq.'run_snowtran') then
            ipar_count = ipar_count + 1
            cpar_name(ipar_count) = c_param(1:i_param_chars)
            call char2real(run_snowtran,i_value_chars,c_value,&
     &        c_param(1:i_param_chars))
            if (run_snowtran.ne.0.0 .and. run_snowtran.ne.1.0) then
              print *,'run_snowtran not 0.0 or 1.0'
              stop
            endif
          endif

          if (c_param(1:i_param_chars).eq.'irun_data_assim') then
            ipar_count = ipar_count + 1
            cpar_name(ipar_count) = c_param(1:i_param_chars)
            call char2int(irun_data_assim,i_value_chars,c_value,&
     &        c_param(1:i_param_chars))
            if (irun_data_assim.ne.0 .and. irun_data_assim.ne.1) then
              print *,'irun_data_assim not 0 or 1'
              stop
            endif
          endif

          if (c_param(1:i_param_chars).eq.'ihrestart_flag') then
            ipar_count = ipar_count + 1
            cpar_name(ipar_count) = c_param(1:i_param_chars)
            call char2int(ihrestart_flag,i_value_chars,c_value,&
     &        c_param(1:i_param_chars))
            if (ihrestart_flag.lt.-2 .or. ihrestart_flag.gt.300000) then
              print *,'ihrestart_flag out of range'
              stop
            endif
          endif

          if (c_param(1:i_param_chars).eq.'i_dataassim_loop') then
            ipar_count = ipar_count + 1
            cpar_name(ipar_count) = c_param(1:i_param_chars)
            call char2int(i_dataassim_loop,i_value_chars,c_value,&
     &        c_param(1:i_param_chars))
            if (i_dataassim_loop.lt.-366 .or. i_dataassim_loop.gt.1)&
     &        then
              print *,'i_dataassim_loop out of range'
              stop
            endif
          endif

          if (c_param(1:i_param_chars).eq.'ihrestart_inc') then
            ipar_count = ipar_count + 1
            cpar_name(ipar_count) = c_param(1:i_param_chars)
            call char2int(ihrestart_inc,i_value_chars,c_value,&
     &        c_param(1:i_param_chars))
            if (ihrestart_inc.lt.0 .or. ihrestart_inc.gt.8784) then
              print *,'ihrestart_inc out of range'
              stop
            endif
          endif

! MICROMET MODEL SETUP.
          if (c_param(1:i_param_chars).eq.'i_tair_flag') then
            ipar_count = ipar_count + 1
            cpar_name(ipar_count) = c_param(1:i_param_chars)
            call char2int(i_tair_flag,i_value_chars,c_value,&
     &        c_param(1:i_param_chars))
            if (i_tair_flag.ne.0 .and. i_tair_flag.ne.1) then
              print *,'i_tair_flag not 0 or 1'
              stop
            endif
          endif

          if (c_param(1:i_param_chars).eq.'i_rh_flag') then
            ipar_count = ipar_count + 1
            cpar_name(ipar_count) = c_param(1:i_param_chars)
            call char2int(i_rh_flag,i_value_chars,c_value,&
     &        c_param(1:i_param_chars))
            if (i_rh_flag.ne.0 .and. i_rh_flag.ne.1) then
              print *,'i_rh_flag not 0 or 1'
              stop
            endif
          endif

          if (c_param(1:i_param_chars).eq.'i_wind_flag') then
            ipar_count = ipar_count + 1
            cpar_name(ipar_count) = c_param(1:i_param_chars)
            call char2int(i_wind_flag,i_value_chars,c_value,&
     &        c_param(1:i_param_chars))
            if (i_wind_flag.ne.-1 .and. i_wind_flag.ne.0 .and. &
     &        i_wind_flag.ne.1) then
              print *,'i_wind_flag not -1, 0, or 1'
              stop
            endif
          endif

          if (c_param(1:i_param_chars).eq.'i_solar_flag') then
            ipar_count = ipar_count + 1
            cpar_name(ipar_count) = c_param(1:i_param_chars)
            call char2int(i_solar_flag,i_value_chars,c_value,&
     &        c_param(1:i_param_chars))
            if (i_solar_flag.ne.0 .and. i_solar_flag.ne.1) then
              print *,'i_solar_flag not 0 or 1'
              stop
            endif
          endif

          if (c_param(1:i_param_chars).eq.'i_longwave_flag') then
            ipar_count = ipar_count + 1
            cpar_name(ipar_count) = c_param(1:i_param_chars)
            call char2int(i_longwave_flag,i_value_chars,c_value,&
     &        c_param(1:i_param_chars))
            if (i_longwave_flag.ne.0 .and. i_longwave_flag.ne.1) then
              print *,'i_longwave_flag not 0 or 1'
              stop
            endif
          endif

          if (c_param(1:i_param_chars).eq.'i_prec_flag') then
            ipar_count = ipar_count + 1
            cpar_name(ipar_count) = c_param(1:i_param_chars)
            call char2int(i_prec_flag,i_value_chars,c_value,&
     &        c_param(1:i_param_chars))
            if (i_prec_flag.ne.0 .and. i_prec_flag.ne.1) then
              print *,'i_prec_flag not 0 or 1'
              stop
            endif
          endif

          if (c_param(1:i_param_chars).eq.'ifill') then
            ipar_count = ipar_count + 1
            cpar_name(ipar_count) = c_param(1:i_param_chars)
            call char2int(ifill,i_value_chars,c_value,&
     &        c_param(1:i_param_chars))
            if (ifill.ne.0 .and. ifill.ne.1) then
              print *,'ifill not 0 or 1'
              stop
            endif
          endif

          if (c_param(1:i_param_chars).eq.'iobsint') then
            ipar_count = ipar_count + 1
            cpar_name(ipar_count) = c_param(1:i_param_chars)
            call char2int(iobsint,i_value_chars,c_value,&
     &        c_param(1:i_param_chars))
            if (iobsint.ne.0 .and. iobsint.ne.1) then
              print *,'iobsint not 0 or 1'
              stop
            endif
          endif

          if (c_param(1:i_param_chars).eq.'dn') then
            ipar_count = ipar_count + 1
            cpar_name(ipar_count) = c_param(1:i_param_chars)
            call char2real(dn,i_value_chars,c_value,&
     &        c_param(1:i_param_chars))
            if (dn.lt.1.0 .or. dn.gt.10000.0) then
              print *,'dn out of range'
              stop
            endif
          endif

          if (c_param(1:i_param_chars).eq.'barnes_lg_domain') then
            ipar_count = ipar_count + 1
            cpar_name(ipar_count) = c_param(1:i_param_chars)
            call char2real(barnes_lg_domain,i_value_chars,c_value,&
     &        c_param(1:i_param_chars))
            if (barnes_lg_domain.ne.0.0 .and. barnes_lg_domain.ne.1.0&
     &        .and. barnes_lg_domain.ne.2.0) then
              print *,'barnes_lg_domain not 0.0 or 1.0 or 2.0'
              stop
            endif
          endif

          if (c_param(1:i_param_chars).eq.'n_stns_used') then
            ipar_count = ipar_count + 1
            cpar_name(ipar_count) = c_param(1:i_param_chars)
            call char2int(n_stns_used,i_value_chars,c_value,&
     &        c_param(1:i_param_chars))
            if (n_stns_used.lt.1 .or. n_stns_used.gt.9) then
              print *,'n_stns_used out of range'
              stop
            endif
          endif

          if (c_param(1:i_param_chars).eq.'snowmodel_line_flag') then
            ipar_count = ipar_count + 1
            cpar_name(ipar_count) = c_param(1:i_param_chars)
            call char2real(snowmodel_line_flag,i_value_chars,c_value,&
     &        c_param(1:i_param_chars))
            if (snowmodel_line_flag.ne.0.0 .and.&
     &        snowmodel_line_flag.ne.1.0)&
     &        then
              print *,'snowmodel_line_flag not 0.0 or 1.0'
              stop
            endif
          endif

          if (c_param(1:i_param_chars).eq.'check_met_data') then
            ipar_count = ipar_count + 1
            cpar_name(ipar_count) = c_param(1:i_param_chars)
            call char2real(check_met_data,i_value_chars,c_value,&
     &        c_param(1:i_param_chars))
            if (check_met_data.ne.0.0 .and. check_met_data.ne.1.0)&
     &        then
              print *,'check_met_data not 0.0 or 1.0'
              stop
            endif
          endif

          if (c_param(1:i_param_chars).eq.'curve_len_scale') then
            ipar_count = ipar_count + 1
            cpar_name(ipar_count) = c_param(1:i_param_chars)
            call char2real(curve_len_scale,i_value_chars,c_value,&
     &        c_param(1:i_param_chars))
            if (curve_len_scale.le.0.0 .or. curve_len_scale.gt.5000.0)&
     &        then
              print *,'curve_len_scale out of range'
              stop
            endif
          endif

          if (c_param(1:i_param_chars).eq.'slopewt') then
            ipar_count = ipar_count + 1
            cpar_name(ipar_count) = c_param(1:i_param_chars)
            call char2real(slopewt,i_value_chars,c_value,&
     &        c_param(1:i_param_chars))
            if (slopewt.lt.0.0 .or. slopewt.gt.1.0) then
              print *,'slopewt out of range'
              stop
            endif
          endif

          if (c_param(1:i_param_chars).eq.'curvewt') then
            ipar_count = ipar_count + 1
            cpar_name(ipar_count) = c_param(1:i_param_chars)
            call char2real(curvewt,i_value_chars,c_value,&
     &        c_param(1:i_param_chars))
            if (curvewt.lt.0.0 .or. curvewt.gt.1.0) then
              print *,'curvewt out of range'
              stop
            endif
          endif

          if (c_param(1:i_param_chars).eq.'curve_lg_scale_flag') then
            ipar_count = ipar_count + 1
            cpar_name(ipar_count) = c_param(1:i_param_chars)
            call char2real(curve_lg_scale_flag,i_value_chars,c_value,&
     &        c_param(1:i_param_chars))
            if (curve_lg_scale_flag.ne.0.0 .and. &
     &        curve_lg_scale_flag.ne.1.0) then
              print *,'curve_lg_scale_flag not 0.0 or 1.0'
              stop
            endif
          endif

          if (c_param(1:i_param_chars).eq.'windspd_min') then
            ipar_count = ipar_count + 1
            cpar_name(ipar_count) = c_param(1:i_param_chars)
            call char2real(windspd_min,i_value_chars,c_value,&
     &        c_param(1:i_param_chars))
            if (windspd_min.lt.0.1 .or. windspd_min.gt.2.0) then
              print *,'windspd_min out of range'
              stop
            endif
          endif

          if (c_param(1:i_param_chars).eq.'lapse_rate_user_flag') then
            ipar_count = ipar_count + 1
            cpar_name(ipar_count) = c_param(1:i_param_chars)
            call char2int(lapse_rate_user_flag,i_value_chars,c_value,&
     &        c_param(1:i_param_chars))
            if (lapse_rate_user_flag.ne.0 .and. &
     &        lapse_rate_user_flag.ne.1) then
              print *,'lapse_rate_user_flag not 0 or 1'
              stop
            endif
          endif

          if (c_param(1:i_param_chars).eq. &
     &      'iprecip_lapse_rate_user_flag') then
            ipar_count = ipar_count + 1
            cpar_name(ipar_count) = c_param(1:i_param_chars)
            call char2int(iprecip_lapse_rate_user_flag,i_value_chars, &
     &        c_value,c_param(1:i_param_chars))
            if (iprecip_lapse_rate_user_flag.ne.0 .and. &
     &        iprecip_lapse_rate_user_flag.ne.1) then
              print *,'iprecip_lapse_rate_user_flag not 0 or 1'
              stop
            endif
          endif

          if (c_param(1:i_param_chars).eq.'iprecip_scheme') then
            ipar_count = ipar_count + 1
            cpar_name(ipar_count) = c_param(1:i_param_chars)
            call char2int(iprecip_scheme,i_value_chars, &
     &        c_value,c_param(1:i_param_chars))
            if (iprecip_scheme.ne.1 .and. iprecip_scheme.ne.2) then
              print *,'iprecip_scheme not 1 or 2'
              stop
            endif
          endif

          if (c_param(1:i_param_chars).eq.'snowfall_frac') then
            ipar_count = ipar_count + 1
            cpar_name(ipar_count) = c_param(1:i_param_chars)
            call char2real(snowfall_frac,i_value_chars,c_value, &
     &        c_param(1:i_param_chars))
            if (snowfall_frac.ne.1.0 .and. snowfall_frac.ne.2.0 &
     &        .and. snowfall_frac.ne.3.0) then
              print *,'snowfall_frac must be = 1.0, 2.0, or 3.0'
              stop
            endif
          endif

          if (c_param(1:i_param_chars).eq.'wind_lapse_rate') then
            ipar_count = ipar_count + 1
            cpar_name(ipar_count) = c_param(1:i_param_chars)
            call char2real(wind_lapse_rate,i_value_chars,c_value, &
     &        c_param(1:i_param_chars))
            if (wind_lapse_rate.lt.0.0) then
              print *,'wind_lapse_rate must be >= 0.0'
              stop
            endif
          endif

          if (c_param(1:i_param_chars).eq.'calc_subcanopy_met') then
            ipar_count = ipar_count + 1
            cpar_name(ipar_count) = c_param(1:i_param_chars)
            call char2real(calc_subcanopy_met,i_value_chars,c_value, &
     &        c_param(1:i_param_chars))
            if (calc_subcanopy_met.ne.0.0 .and. &
     &        calc_subcanopy_met.ne.1.0) then
              print *,'calc_subcanopy_met not 0.0 or 1.0'
              stop
            endif
          endif

          if (c_param(1:i_param_chars).eq.'gap_frac') then
            ipar_count = ipar_count + 1
            cpar_name(ipar_count) = c_param(1:i_param_chars)
            call char2real(gap_frac,i_value_chars,c_value, &
     &        c_param(1:i_param_chars))
            if (gap_frac.lt.0.0 .or. gap_frac.gt.1.0) then
              print *,'gap_frac out of range'
              stop
            endif
          endif

          if (c_param(1:i_param_chars).eq.'cloud_frac_factor') then
            ipar_count = ipar_count + 1
            cpar_name(ipar_count) = c_param(1:i_param_chars)
            call char2real(cloud_frac_factor,i_value_chars,c_value, &
     &        c_param(1:i_param_chars))
            if (cloud_frac_factor.lt.0.0 .or. cloud_frac_factor.gt.1.0) &
     &        then
              print *,'cloud_frac_factor out of range'
              stop
            endif
          endif

          if (c_param(1:i_param_chars).eq.'use_shortwave_obs') then
            ipar_count = ipar_count + 1
            cpar_name(ipar_count) = c_param(1:i_param_chars)
            call char2real(use_shortwave_obs,i_value_chars,c_value,&
     &        c_param(1:i_param_chars))
            if (use_shortwave_obs.ne.0.0 .and. use_shortwave_obs.ne.1.0) &
     &        then
              print *,'use_shortwave_obs not 0.0 or 1.0'
              stop
            endif
          endif

          if (c_param(1:i_param_chars).eq.'use_longwave_obs') then
            ipar_count = ipar_count + 1
            cpar_name(ipar_count) = c_param(1:i_param_chars)
            call char2real(use_longwave_obs,i_value_chars,c_value,&
     &        c_param(1:i_param_chars))
            if (use_longwave_obs.ne.0.0 .and. use_longwave_obs.ne.1.0) &
     &        then
              print *,'use_longwave_obs not 0.0 or 1.0'
              stop
            endif
          endif

          if (c_param(1:i_param_chars).eq.'use_sfc_pressure_obs') then
            ipar_count = ipar_count + 1
            cpar_name(ipar_count) = c_param(1:i_param_chars)
            call char2real(use_sfc_pressure_obs,i_value_chars,c_value,&
     &        c_param(1:i_param_chars))
            if (use_sfc_pressure_obs.ne.0.0 .and. &
     &        use_sfc_pressure_obs.ne.1.0) then
              print *,'use_sfc_pressure_obs not 0.0 or 1.0'
              stop
            endif
          endif

          if (c_param(1:i_param_chars).eq.'cf_precip_flag') then
            ipar_count = ipar_count + 1
            cpar_name(ipar_count) = c_param(1:i_param_chars)
            call char2real(cf_precip_flag,i_value_chars,c_value,&
     &        c_param(1:i_param_chars))
            if (cf_precip_flag.ne.0.0 .and. cf_precip_flag.ne.1.0 .and. &
     &        cf_precip_flag.ne.2.0 .and. cf_precip_flag.ne.3.0) then
              print *,'cf_precip_flag not 0.0, 1.0, 2.0, or 3.0'
              stop
            endif
          endif

! SNOWTRAN-3D MODEL SETUP.
          if (c_param(1:i_param_chars).eq.'Utau_t_flag') then
            ipar_count = ipar_count + 1
            cpar_name(ipar_count) = c_param(1:i_param_chars)
            call char2real(Utau_t_flag,i_value_chars,c_value,&
     &        c_param(1:i_param_chars))
            if (Utau_t_flag.ne.0.0 .and. Utau_t_flag.ne.1.0) then
              print *,'Utau_t_flag not 0.0 or 1.0'
              stop
            endif
          endif

          if (c_param(1:i_param_chars).eq.'Utau_t_const') then
            ipar_count = ipar_count + 1
            cpar_name(ipar_count) = c_param(1:i_param_chars)
            call char2real(Utau_t_const,i_value_chars,c_value,&
     &        c_param(1:i_param_chars))
            if (Utau_t_const.le.0.0 .and. Utau_t_const.gt.1.0) then
              print *,'Utau_t_const out of range'
              stop
            endif
          endif

          if (c_param(1:i_param_chars).eq.'subgrid_flag') then
            ipar_count = ipar_count + 1
            cpar_name(ipar_count) = c_param(1:i_param_chars)
            call char2real(subgrid_flag,i_value_chars,c_value,&
     &        c_param(1:i_param_chars))
            if (subgrid_flag.ne.0.0 .and. subgrid_flag.ne.1.0) then
              print *,'subgrid_flag not 0.0 or 1.0'
              stop
            endif
          endif

          if (c_param(1:i_param_chars).eq.'tabler_dir') then
            ipar_count = ipar_count + 1
            cpar_name(ipar_count) = c_param(1:i_param_chars)
            call char2real(tabler_dir,i_value_chars,c_value,&
     &        c_param(1:i_param_chars))
            if (tabler_dir.lt.-360.0 .or. tabler_dir.gt.360.0) then
              print *,'tabler_dir out of range'
              stop
            endif
          endif

          if (c_param(1:i_param_chars).eq.'slope_adjust') then
            ipar_count = ipar_count + 1
            cpar_name(ipar_count) = c_param(1:i_param_chars)
            call char2real(slope_adjust,i_value_chars,c_value,&
     &        c_param(1:i_param_chars))
            if (slope_adjust.lt.0.0 .or. slope_adjust.gt.3.0) then
              print *,'slope_adjust out of range'
              stop
            endif
          endif

          if (c_param(1:i_param_chars).eq.'twolayer_flag') then
            ipar_count = ipar_count + 1
            cpar_name(ipar_count) = c_param(1:i_param_chars)
            call char2real(twolayer_flag,i_value_chars,c_value,&
     &        c_param(1:i_param_chars))
            if (twolayer_flag.ne.0.0 .and. twolayer_flag.ne.1.0) then
              print *,'twolayer_flag not 0.0 or 1.0'
              stop
            endif
          endif

          if (c_param(1:i_param_chars).eq.'bc_flag') then
            ipar_count = ipar_count + 1
            cpar_name(ipar_count) = c_param(1:i_param_chars)
            call char2real(bc_flag,i_value_chars,c_value,&
     &        c_param(1:i_param_chars))
            if (bc_flag.ne.0.0 .and. bc_flag.ne.1.0) then
              print *,'bc_flag not 0.0 or 1.0'
              stop
            endif
          endif

          if (c_param(1:i_param_chars).eq.'ht_windobs') then
            ipar_count = ipar_count + 1
            cpar_name(ipar_count) = c_param(1:i_param_chars)
            call char2real(ht_windobs,i_value_chars,c_value,&
     &        c_param(1:i_param_chars))
            if (ht_windobs.le.0.0 .or. ht_windobs.gt.20.0) then
              print *,'ht_windobs out of range'
              stop
            endif
          endif

          if (c_param(1:i_param_chars).eq.'ht_rhobs') then
            ipar_count = ipar_count + 1
            cpar_name(ipar_count) = c_param(1:i_param_chars)
            call char2real(ht_rhobs,i_value_chars,c_value,&
     &        c_param(1:i_param_chars))
            if (ht_rhobs.le.0.0 .or. ht_rhobs.gt.20.0) then
              print *,'ht_rhobs out of range'
              stop
            endif
          endif

          if (c_param(1:i_param_chars).eq.'ro_snow') then
            ipar_count = ipar_count + 1
            cpar_name(ipar_count) = c_param(1:i_param_chars)
            call char2real(ro_snow,i_value_chars,c_value,&
     &        c_param(1:i_param_chars))
            if (ro_snow.le.100.0 .or. ro_snow.gt.550.0) then
              print *,'ro_snow out of range'
              stop
            endif
          endif

          if (c_param(1:i_param_chars).eq.'snow_d_init_const') then
            ipar_count = ipar_count + 1
            cpar_name(ipar_count) = c_param(1:i_param_chars)
            call char2real(snow_d_init_const,i_value_chars,c_value,&
     &        c_param(1:i_param_chars))
            if (snow_d_init_const.lt.0.0 .or. &
     &        snow_d_init_const.gt.5.0) then
              print *,'snow_d_init_const out of range'
              stop
            endif
          endif

          if (c_param(1:i_param_chars).eq.'topoflag') then
            ipar_count = ipar_count + 1
            cpar_name(ipar_count) = c_param(1:i_param_chars)
            call char2real(topoflag,i_value_chars,c_value,&
     &        c_param(1:i_param_chars))
            if (topoflag.ne.0.0 .and. topoflag.ne.1.0) then
              print *,'topoflag not 0.0 or 1.0'
              stop
            endif
          endif

! ENBAL MODEL SETUP.
          if (c_param(1:i_param_chars).eq.'icond_flag') then
            ipar_count = ipar_count + 1
            cpar_name(ipar_count) = c_param(1:i_param_chars)
            call char2int(icond_flag,i_value_chars,c_value,&
     &        c_param(1:i_param_chars))
            if (icond_flag.ne.0 .and. icond_flag.ne.1) then
              print *,'icond_flag not 0 or 1'
              stop
            endif
          endif

          if (c_param(1:i_param_chars).eq.'albedo_snow_forest') then
            ipar_count = ipar_count + 1
            cpar_name(ipar_count) = c_param(1:i_param_chars)
            call char2real(albedo_snow_forest,i_value_chars,c_value,&
     &        c_param(1:i_param_chars))
            if (albedo_snow_forest.lt.0.3 .or. &
     &        albedo_snow_forest.gt.0.8) then
              print *,'albedo_snow_forest out of range'
              stop
            endif
          endif

          if (c_param(1:i_param_chars).eq.'albedo_snow_clearing') then
            ipar_count = ipar_count + 1
            cpar_name(ipar_count) = c_param(1:i_param_chars)
            call char2real(albedo_snow_clearing,i_value_chars,c_value,&
     &        c_param(1:i_param_chars))
            if (albedo_snow_clearing.lt.0.4 .or. &
     &        albedo_snow_clearing.gt.0.8) then
              print *,'albedo_snow_clearing out of range'
              stop
            endif
          endif

          if (c_param(1:i_param_chars).eq.'albedo_glacier') then
            ipar_count = ipar_count + 1
            cpar_name(ipar_count) = c_param(1:i_param_chars)
            call char2real(albedo_glacier,i_value_chars,c_value,&
     &        c_param(1:i_param_chars))
            if (albedo_glacier.lt.0.3 .or. albedo_glacier.gt.0.8) then
              print *,'albedo_glacier out of range'
              stop
            endif
          endif

! SNOWPACK MODEL SETUP.
          if (c_param(1:i_param_chars).eq.'sfc_sublim_flag') then
            ipar_count = ipar_count + 1
            cpar_name(ipar_count) = c_param(1:i_param_chars)
            call char2real(sfc_sublim_flag,i_value_chars,c_value,&
     &        c_param(1:i_param_chars))
            if (sfc_sublim_flag.ne.0.0 .and. sfc_sublim_flag.ne.1.0)&
     &        then
              print *,'sfc_sublim_flag not 0.0 or 1.0'
              stop
            endif
          endif

          if (c_param(1:i_param_chars).eq.'multilayer_snowpack') then
            ipar_count = ipar_count + 1
            cpar_name(ipar_count) = c_param(1:i_param_chars)
            call char2int(multilayer_snowpack,i_value_chars,c_value,&
     &        c_param(1:i_param_chars))
            if (multilayer_snowpack.ne.0 .and. multilayer_snowpack.ne.1)&
     &        then
              print *,'multilayer_snowpack not 0 or 1'
              stop
            endif
          endif

          if (c_param(1:i_param_chars).eq.'tsls_threshold') then
            ipar_count = ipar_count + 1
            cpar_name(ipar_count) = c_param(1:i_param_chars)
            call char2real(tsls_threshold,i_value_chars,c_value,&
     &        c_param(1:i_param_chars))
            if (tsls_threshold.lt.1.0 .or. tsls_threshold.gt.8760.0) &
     &        then
              print *,'tsls_threshold out of range'
              stop
            endif
          endif

          if (c_param(1:i_param_chars).eq.'max_layers') then
            ipar_count = ipar_count + 1
            cpar_name(ipar_count) = c_param(1:i_param_chars)
            call char2int(max_layers,i_value_chars,c_value,&
     &        c_param(1:i_param_chars))
            if (max_layers.lt.1 .or. max_layers.gt.100) then
              print *,'max_layers out of range'
              stop
            endif
          endif

          if (c_param(1:i_param_chars).eq.'dz_snow_min') then
            ipar_count = ipar_count + 1
            cpar_name(ipar_count) = c_param(1:i_param_chars)
            call char2real(dz_snow_min,i_value_chars,c_value,&
     &        c_param(1:i_param_chars))
            if (dz_snow_min.lt.0.0 .or. dz_snow_min.gt.5.0) then
              print *,'dz_snow_min out of range'
              stop
            endif
          endif

          if (c_param(1:i_param_chars).eq.'izero_snow_date') then
            ipar_count = ipar_count + 1
            cpar_name(ipar_count) = c_param(1:i_param_chars)
            call char2int(izero_snow_date,i_value_chars,c_value,&
     &        c_param(1:i_param_chars))
            if (izero_snow_date.lt.10100 .or. &
     &        izero_snow_date.gt.999999) then
              print *,'izero_snow_date out of range'
              stop
            endif
          endif

! SEAICE MODEL SETUP.
          if (c_param(1:i_param_chars).eq.'seaice_run') then
            ipar_count = ipar_count + 1
            cpar_name(ipar_count) = c_param(1:i_param_chars)
            call char2real(seaice_run,i_value_chars,c_value,&
     &        c_param(1:i_param_chars))
            if (seaice_run.ne.0.0 .and. seaice_run.ne.1.0 .and. &
     &        seaice_run.ne.2.0 .and. seaice_run.ne.3.0 .and. &
     &        seaice_run.ne.4.0) then
              print *,'seaice_run not 0.0, 1.0, 2.0, 3.0, or 4.0'
              stop
            endif
          endif

! PRINTING OUTPUT FILES.
          if (c_param(1:i_param_chars).eq.'print_micromet') then
            ipar_count = ipar_count + 1
            cpar_name(ipar_count) = c_param(1:i_param_chars)
            call char2real(print_micromet,i_value_chars,c_value,&
     &        c_param(1:i_param_chars))
            if (print_micromet.ne.0.0 .and. print_micromet.ne.1.0) then
              print *,'print_micromet not 0.0 or 1.0'
              stop
            endif
          endif

          if (c_param(1:i_param_chars).eq.'micromet_output_fname') then
            ipar_count = ipar_count + 1
            cpar_name(ipar_count) = c_param(1:i_param_chars)
            call char2char(micromet_output_fname,c_value,&
     &        i_value_chars,c_param(1:i_param_chars))
          endif

          if (c_param(1:i_param_chars).eq.'print_snowtran') then
            ipar_count = ipar_count + 1
            cpar_name(ipar_count) = c_param(1:i_param_chars)
            call char2real(print_snowtran,i_value_chars,c_value,&
     &        c_param(1:i_param_chars))
            if (print_snowtran.ne.0.0 .and. print_snowtran.ne.1.0) then
              print *,'print_snowtran not 0.0 or 1.0'
              stop
            endif
          endif

          if (c_param(1:i_param_chars).eq.'snowtran_output_fname') then
            ipar_count = ipar_count + 1
            cpar_name(ipar_count) = c_param(1:i_param_chars)
            call char2char(snowtran_output_fname,c_value,&
     &        i_value_chars,c_param(1:i_param_chars))
          endif

          if (c_param(1:i_param_chars).eq.'Tabler_1_flag') then
            ipar_count = ipar_count + 1
            cpar_name(ipar_count) = c_param(1:i_param_chars)
            call char2real(Tabler_1_flag,i_value_chars,c_value,&
     &        c_param(1:i_param_chars))
            if (Tabler_1_flag.ne.0.0 .and. Tabler_1_flag.ne.1.0) then
              print *,'Tabler_1_flag not 0.0 or 1.0'
              stop
            endif
          endif

          if (c_param(1:i_param_chars).eq.'tabler_sfc_path_name') then
            ipar_count = ipar_count + 1
            cpar_name(ipar_count) = c_param(1:i_param_chars)
            call char2char(tabler_sfc_path_name,c_value,&
     &        i_value_chars,c_param(1:i_param_chars))
          endif

          if (c_param(1:i_param_chars).eq.'Tabler_2_flag') then
            ipar_count = ipar_count + 1
            cpar_name(ipar_count) = c_param(1:i_param_chars)
            call char2real(Tabler_2_flag,i_value_chars,c_value,&
     &        c_param(1:i_param_chars))
            if (Tabler_2_flag.ne.0.0 .and. Tabler_2_flag.ne.1.0) then
              print *,'Tabler_2_flag not 0.0 or 1.0'
              stop
            endif
          endif

          if (c_param(1:i_param_chars).eq.'print_enbal') then
            ipar_count = ipar_count + 1
            cpar_name(ipar_count) = c_param(1:i_param_chars)
            call char2real(print_enbal,i_value_chars,c_value,&
     &        c_param(1:i_param_chars))
            if (print_enbal.ne.0.0 .and. print_enbal.ne.1.0) then
              print *,'print_enbal not 0.0 or 1.0'
              stop
            endif
          endif

          if (c_param(1:i_param_chars).eq.'enbal_output_fname') then
            ipar_count = ipar_count + 1
            cpar_name(ipar_count) = c_param(1:i_param_chars)
            call char2char(enbal_output_fname,c_value,&
     &        i_value_chars,c_param(1:i_param_chars))
          endif

          if (c_param(1:i_param_chars).eq.'print_snowpack') then
            ipar_count = ipar_count + 1
            cpar_name(ipar_count) = c_param(1:i_param_chars)
            call char2real(print_snowpack,i_value_chars,c_value,&
     &        c_param(1:i_param_chars))
            if (print_snowpack.ne.0.0 .and. print_snowpack.ne.1.0) then
              print *,'print_snowpack not 0.0 or 1.0'
              stop
            endif
          endif

          if (c_param(1:i_param_chars).eq.'snowpack_output_fname') then
            ipar_count = ipar_count + 1
            cpar_name(ipar_count) = c_param(1:i_param_chars)
            call char2char(snowpack_output_fname,c_value,&
     &        i_value_chars,c_param(1:i_param_chars))
          endif

          if (c_param(1:i_param_chars).eq.'print_multilayer') then
            ipar_count = ipar_count + 1
            cpar_name(ipar_count) = c_param(1:i_param_chars)
            call char2real(print_multilayer,i_value_chars,c_value,&
     &        c_param(1:i_param_chars))
            if (print_multilayer.ne.0.0 .and. print_multilayer.ne.1.0&
     &        .and. print_multilayer.ne.2.0) &
     &        then
              print *,'print_multilayer not 0.0 or 1.0 or 2.0'
              stop
            endif
          endif

          if (c_param(1:i_param_chars).eq. &
     &      'multilayer_output_fname') then
            ipar_count = ipar_count + 1
            cpar_name(ipar_count) = c_param(1:i_param_chars)
            call char2char(multilayer_output_fname,c_value,&
     &        i_value_chars,c_param(1:i_param_chars))
          endif

          if (c_param(1:i_param_chars).eq.'print_user') then
            ipar_count = ipar_count + 1
            cpar_name(ipar_count) = c_param(1:i_param_chars)
            call char2real(print_user,i_value_chars,c_value,&
     &        c_param(1:i_param_chars))
            if (print_user.ne.0.0 .and. print_user.ne.1.0) then
              print *,'print_user not 0.0 or 1.0'
              stop
            endif
          endif

          if (c_param(1:i_param_chars).eq.'print_inc') then
            ipar_count = ipar_count + 1
            cpar_name(ipar_count) = c_param(1:i_param_chars)
            call char2real(print_inc,i_value_chars,c_value,&
     &        c_param(1:i_param_chars))
            if (print_inc.ne.1.0 .and. print_inc.ne.8.0 .and.&
     &       print_inc.ne.24.0) then
              print *,'print_inc not 1.0 or 8.0 or 24.0.  You are'
              print *,'using some number that is not typical and'
              print *,'should make certain you and the code are' 
              print *,'doing exactly what you think it is and want'
              print *,'to do.'
              stop
            endif
          endif

          if (c_param(1:i_param_chars).eq.'output_path_wo_assim') then
            ipar_count = ipar_count + 1
            cpar_name(ipar_count) = c_param(1:i_param_chars)
            call char2char(output_path_wo_assim,c_value,&
     &        i_value_chars,c_param(1:i_param_chars))
          endif

          if (c_param(1:i_param_chars).eq.'output_path_wi_assim') then
            ipar_count = ipar_count + 1
            cpar_name(ipar_count) = c_param(1:i_param_chars)
            call char2char(output_path_wi_assim,c_value,&
     &        i_value_chars,c_param(1:i_param_chars))
          endif

! THESE ARE THE 30 PRINT_VAR_## FLAGS.
          if (c_param(1:i_param_chars).eq.'print_var_01') then
            ipar_count = ipar_count + 1
            cpar_name(ipar_count) = c_param(1:i_param_chars)
            call char2char(print_var_01,c_value, &
     &        i_value_chars,c_param(1:i_param_chars))
            if (print_var_01.ne.'y' .and. print_var_01.ne.'n') then
              print *,'print_var_01 not y or n'
              stop
            endif
          endif

          if (c_param(1:i_param_chars).eq.'print_var_02') then
            ipar_count = ipar_count + 1
            cpar_name(ipar_count) = c_param(1:i_param_chars)
            call char2char(print_var_02,c_value,&
     &        i_value_chars,c_param(1:i_param_chars))
            if (print_var_02.ne.'y' .and. print_var_02.ne.'n') then
              print *,'print_var_02 not y or n'
              stop
            endif
          endif

          if (c_param(1:i_param_chars).eq.'print_var_03') then
            ipar_count = ipar_count + 1
            cpar_name(ipar_count) = c_param(1:i_param_chars)
            call char2char(print_var_03,c_value,&
     &        i_value_chars,c_param(1:i_param_chars))
            if (print_var_03.ne.'y' .and. print_var_03.ne.'n') then
              print *,'print_var_03 not y or n'
              stop
            endif
          endif

          if (c_param(1:i_param_chars).eq.'print_var_04') then
            ipar_count = ipar_count + 1
            cpar_name(ipar_count) = c_param(1:i_param_chars)
            call char2char(print_var_04,c_value,&
     &        i_value_chars,c_param(1:i_param_chars))
            if (print_var_04.ne.'y' .and. print_var_04.ne.'n') then
              print *,'print_var_04 not y or n'
              stop
            endif
          endif

          if (c_param(1:i_param_chars).eq.'print_var_05') then
            ipar_count = ipar_count + 1
            cpar_name(ipar_count) = c_param(1:i_param_chars)
            call char2char(print_var_05,c_value,&
     &        i_value_chars,c_param(1:i_param_chars))
            if (print_var_05.ne.'y' .and. print_var_05.ne.'n') then
              print *,'print_var_05 not y or n'
              stop
            endif
          endif

          if (c_param(1:i_param_chars).eq.'print_var_06') then
            ipar_count = ipar_count + 1
            cpar_name(ipar_count) = c_param(1:i_param_chars)
            call char2char(print_var_06,c_value,&
     &        i_value_chars,c_param(1:i_param_chars))
            if (print_var_06.ne.'y' .and. print_var_06.ne.'n') then
              print *,'print_var_06 not y or n'
              stop
            endif
          endif

          if (c_param(1:i_param_chars).eq.'print_var_07') then
            ipar_count = ipar_count + 1
            cpar_name(ipar_count) = c_param(1:i_param_chars)
            call char2char(print_var_07,c_value,&
     &        i_value_chars,c_param(1:i_param_chars))
            if (print_var_07.ne.'y' .and. print_var_07.ne.'n') then
              print *,'print_var_07 not y or n'
              stop
            endif
          endif

          if (c_param(1:i_param_chars).eq.'print_var_08') then
            ipar_count = ipar_count + 1
            cpar_name(ipar_count) = c_param(1:i_param_chars)
            call char2char(print_var_08,c_value,&
     &        i_value_chars,c_param(1:i_param_chars))
            if (print_var_08.ne.'y' .and. print_var_08.ne.'n') then
              print *,'print_var_08 not y or n'
              stop
            endif
          endif

          if (c_param(1:i_param_chars).eq.'print_var_09') then
            ipar_count = ipar_count + 1
            cpar_name(ipar_count) = c_param(1:i_param_chars)
            call char2char(print_var_09,c_value,&
     &        i_value_chars,c_param(1:i_param_chars))
            if (print_var_09.ne.'y' .and. print_var_09.ne.'n') then
              print *,'print_var_09 not y or n'
              stop
            endif
          endif

          if (c_param(1:i_param_chars).eq.'print_var_10') then
            ipar_count = ipar_count + 1
            cpar_name(ipar_count) = c_param(1:i_param_chars)
            call char2char(print_var_10,c_value,&
     &        i_value_chars,c_param(1:i_param_chars))
            if (print_var_10.ne.'y' .and. print_var_10.ne.'n') then
              print *,'print_var_10 not y or n'
              stop
            endif
          endif

          if (c_param(1:i_param_chars).eq.'print_var_11') then
            ipar_count = ipar_count + 1
            cpar_name(ipar_count) = c_param(1:i_param_chars)
            call char2char(print_var_11,c_value,&
     &        i_value_chars,c_param(1:i_param_chars))
            if (print_var_11.ne.'y' .and. print_var_11.ne.'n') then
              print *,'print_var_11 not y or n'
              stop
            endif
          endif

          if (c_param(1:i_param_chars).eq.'print_var_12') then
            ipar_count = ipar_count + 1
            cpar_name(ipar_count) = c_param(1:i_param_chars)
            call char2char(print_var_12,c_value,&
     &        i_value_chars,c_param(1:i_param_chars))
            if (print_var_12.ne.'y' .and. print_var_12.ne.'n') then
              print *,'print_var_12 not y or n'
              stop
            endif
          endif

          if (c_param(1:i_param_chars).eq.'print_var_13') then
            ipar_count = ipar_count + 1
            cpar_name(ipar_count) = c_param(1:i_param_chars)
            call char2char(print_var_13,c_value,&
     &        i_value_chars,c_param(1:i_param_chars))
            if (print_var_13.ne.'y' .and. print_var_13.ne.'n') then
              print *,'print_var_13 not y or n'
              stop
            endif
          endif

          if (c_param(1:i_param_chars).eq.'print_var_14') then
            ipar_count = ipar_count + 1
            cpar_name(ipar_count) = c_param(1:i_param_chars)
            call char2char(print_var_14,c_value,&
     &        i_value_chars,c_param(1:i_param_chars))
            if (print_var_14.ne.'y' .and. print_var_14.ne.'n') then
              print *,'print_var_14 not y or n'
              stop
            endif
          endif

          if (c_param(1:i_param_chars).eq.'print_var_15') then
            ipar_count = ipar_count + 1
            cpar_name(ipar_count) = c_param(1:i_param_chars)
            call char2char(print_var_15,c_value,&
     &        i_value_chars,c_param(1:i_param_chars))
            if (print_var_15.ne.'y' .and. print_var_15.ne.'n') then
              print *,'print_var_15 not y or n'
              stop
            endif
          endif

          if (c_param(1:i_param_chars).eq.'print_var_16') then
            ipar_count = ipar_count + 1
            cpar_name(ipar_count) = c_param(1:i_param_chars)
            call char2char(print_var_16,c_value,&
     &        i_value_chars,c_param(1:i_param_chars))
            if (print_var_16.ne.'y' .and. print_var_16.ne.'n') then
              print *,'print_var_16 not y or n'
              stop
            endif
          endif

          if (c_param(1:i_param_chars).eq.'print_var_17') then
            ipar_count = ipar_count + 1
            cpar_name(ipar_count) = c_param(1:i_param_chars)
            call char2char(print_var_17,c_value,&
     &        i_value_chars,c_param(1:i_param_chars))
            if (print_var_17.ne.'y' .and. print_var_17.ne.'n') then
              print *,'print_var_17 not y or n'
              stop
            endif
          endif

          if (c_param(1:i_param_chars).eq.'print_var_18') then
            ipar_count = ipar_count + 1
            cpar_name(ipar_count) = c_param(1:i_param_chars)
            call char2char(print_var_18,c_value,&
     &        i_value_chars,c_param(1:i_param_chars))
            if (print_var_18.ne.'y' .and. print_var_18.ne.'n') then
              print *,'print_var_18 not y or n'
              stop
            endif
          endif

          if (c_param(1:i_param_chars).eq.'print_var_19') then
            ipar_count = ipar_count + 1
            cpar_name(ipar_count) = c_param(1:i_param_chars)
            call char2char(print_var_19,c_value,&
     &        i_value_chars,c_param(1:i_param_chars))
            if (print_var_19.ne.'y' .and. print_var_19.ne.'n') then
              print *,'print_var_19 not y or n'
              stop
            endif
          endif

          if (c_param(1:i_param_chars).eq.'print_var_20') then
            ipar_count = ipar_count + 1
            cpar_name(ipar_count) = c_param(1:i_param_chars)
            call char2char(print_var_20,c_value,&
     &        i_value_chars,c_param(1:i_param_chars))
            if (print_var_20.ne.'y' .and. print_var_20.ne.'n') then
              print *,'print_var_20 not y or n'
              stop
            endif
          endif

          if (c_param(1:i_param_chars).eq.'print_var_21') then
            ipar_count = ipar_count + 1
            cpar_name(ipar_count) = c_param(1:i_param_chars)
            call char2char(print_var_21,c_value,&
     &        i_value_chars,c_param(1:i_param_chars))
            if (print_var_21.ne.'y' .and. print_var_21.ne.'n') then
              print *,'print_var_21 not y or n'
              stop
            endif
          endif

          if (c_param(1:i_param_chars).eq.'print_var_22') then
            ipar_count = ipar_count + 1
            cpar_name(ipar_count) = c_param(1:i_param_chars)
            call char2char(print_var_22,c_value,&
     &        i_value_chars,c_param(1:i_param_chars))
            if (print_var_22.ne.'y' .and. print_var_22.ne.'n') then
              print *,'print_var_22 not y or n'
              stop
            endif
          endif

          if (c_param(1:i_param_chars).eq.'print_var_23') then
            ipar_count = ipar_count + 1
            cpar_name(ipar_count) = c_param(1:i_param_chars)
            call char2char(print_var_23,c_value,&
     &        i_value_chars,c_param(1:i_param_chars))
            if (print_var_23.ne.'y' .and. print_var_23.ne.'n') then
              print *,'print_var_23 not y or n'
              stop
            endif
          endif

          if (c_param(1:i_param_chars).eq.'print_var_24') then
            ipar_count = ipar_count + 1
            cpar_name(ipar_count) = c_param(1:i_param_chars)
            call char2char(print_var_24,c_value,&
     &        i_value_chars,c_param(1:i_param_chars))
            if (print_var_24.ne.'y' .and. print_var_24.ne.'n') then
              print *,'print_var_24 not y or n'
              stop
            endif
          endif

          if (c_param(1:i_param_chars).eq.'print_var_25') then
            ipar_count = ipar_count + 1
            cpar_name(ipar_count) = c_param(1:i_param_chars)
            call char2char(print_var_25,c_value,&
     &        i_value_chars,c_param(1:i_param_chars))
            if (print_var_25.ne.'y' .and. print_var_25.ne.'n') then
              print *,'print_var_25 not y or n'
              stop
            endif
          endif

          if (c_param(1:i_param_chars).eq.'print_var_26') then
            ipar_count = ipar_count + 1
            cpar_name(ipar_count) = c_param(1:i_param_chars)
            call char2char(print_var_26,c_value,&
     &        i_value_chars,c_param(1:i_param_chars))
            if (print_var_26.ne.'y' .and. print_var_26.ne.'n') then
              print *,'print_var_26 not y or n'
              stop
            endif
          endif

          if (c_param(1:i_param_chars).eq.'print_var_27') then
            ipar_count = ipar_count + 1
            cpar_name(ipar_count) = c_param(1:i_param_chars)
            call char2char(print_var_27,c_value,&
     &        i_value_chars,c_param(1:i_param_chars))
            if (print_var_27.ne.'y' .and. print_var_27.ne.'n') then
              print *,'print_var_27 not y or n'
              stop
            endif
          endif

          if (c_param(1:i_param_chars).eq.'print_var_28') then
            ipar_count = ipar_count + 1
            cpar_name(ipar_count) = c_param(1:i_param_chars)
            call char2char(print_var_28,c_value,&
     &        i_value_chars,c_param(1:i_param_chars))
            if (print_var_28.ne.'y' .and. print_var_28.ne.'n') then
              print *,'print_var_28 not y or n'
              stop
            endif
          endif

          if (c_param(1:i_param_chars).eq.'print_var_29') then
            ipar_count = ipar_count + 1
            cpar_name(ipar_count) = c_param(1:i_param_chars)
            call char2char(print_var_29,c_value,&
     &        i_value_chars,c_param(1:i_param_chars))
            if (print_var_29.ne.'y' .and. print_var_29.ne.'n') then
              print *,'print_var_29 not y or n'
              stop
            endif
          endif

          if (c_param(1:i_param_chars).eq.'print_var_30') then
            ipar_count = ipar_count + 1
            cpar_name(ipar_count) = c_param(1:i_param_chars)
            call char2char(print_var_30,c_value,&
     &        i_value_chars,c_param(1:i_param_chars))
            if (print_var_30.ne.'y' .and. print_var_30.ne.'n') then
              print *,'print_var_30 not y or n'
              stop
            endif
          endif

! Real example
!         if (c_param(1:i_param_chars).eq.'')
!    &      call char2real(,i_value_chars,c_value,
!    &        c_param(1:i_param_chars))

! Integer example.
!         if (c_param(1:i_param_chars).eq.'nx')
!    &      call char2int(nx,i_value_chars,c_value,
!    &        c_param(1:i_param_chars))

! Character example.
!         if (c_param(1:i_param_chars).eq.'fname_out')
!    &      call char2char(fname_out,c_value,i_value_chars,
!    &        c_param(1:i_param_chars))

        endif
      enddo

  99  continue

! Use the print_var_## parameters to build the print control
!   array "print_var(n_print_vars)" used in output_user.f.
      print_var(1) = print_var_01
      print_var(2) = print_var_02
      print_var(3) = print_var_03
      print_var(4) = print_var_04
      print_var(5) = print_var_05
      print_var(6) = print_var_06
      print_var(7) = print_var_07
      print_var(8) = print_var_08
      print_var(9) = print_var_09
      print_var(10) = print_var_10
      print_var(11) = print_var_11
      print_var(12) = print_var_12
      print_var(13) = print_var_13
      print_var(14) = print_var_14
      print_var(15) = print_var_15
      print_var(16) = print_var_16
      print_var(17) = print_var_17
      print_var(18) = print_var_18
      print_var(19) = print_var_19
      print_var(20) = print_var_20
      print_var(21) = print_var_21
      print_var(22) = print_var_22
      print_var(23) = print_var_23
      print_var(24) = print_var_24
      print_var(25) = print_var_25
      print_var(26) = print_var_26
      print_var(27) = print_var_27
      print_var(28) = print_var_28
      print_var(29) = print_var_29
      print_var(30) = print_var_30

! Check the input-parameter counting array to be sure that all
!   parameters have been defined.
      ipar_flag = npars - ipar_count
      if (ipar_flag.ne.0) then
        print *
        print *
        print *,'THERE ARE MISSING VARIABLES IN THE .PAR FILE.'
        print *
        print *,'  THIS MANY VARIABLES ARE MISSING:',ipar_flag
        print *
        print *,'  THESE ARE THE VARIABLES THAT HAVE BEEN READ IN:'

        do k=1,npars
          print *,'VARIABLE NAME =',k,'   ',cpar_name(k)
        enddo

        print *
        print *,'ALL .PAR VARIABLES MUST BE DEFINED BEFORE'
        print *,'  YOU CAN CONTINUE.'
        print *

        stop
      endif

! Put some space at the end of the parameter printouts.
      print *

      return
      end subroutine READPARAM_CODE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine char2char(cvar,c_value,i_value_chars,&
     &  c_param)

      implicit none

      character*(*) c_value,c_param,cvar
      integer i_value_chars

      cvar = c_value(1:i_value_chars)

!KRA
!      print *,c_param,' = ',cvar(1:i_value_chars)

      return
      end subroutine char2char

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine char2real(rvar,i_value_chars,c_value,&
     &  c_param)

      implicit none

      character*(*) c_value,c_param
      integer i_value_chars
      real rvar
      character*8 form

! Read an real value (rvar) from the character string (c_value).
      write (form,90) i_value_chars
   90 format ('(f',i2,'.0)')
      read (c_value,form) rvar

!KRA
!      print *,c_param,' =',rvar

      return
      end subroutine char2real

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine char2double(dvar,i_value_chars,c_value,&
     &  c_param)

      implicit none

      character*(*) c_value,c_param
      integer i_value_chars
      double precision dvar
      character*8 form

! Read an double precision value (dvar) from the character
!   string (c_value).
      write (form,90) i_value_chars
   90 format ('(f',i2,'.0)')
      read (c_value,form) dvar

!KRA
!      print *,c_param,' =',dvar

      return
      end subroutine char2double

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine char2int(ivar,i_value_chars,c_value,&
     &  c_param)

      implicit none

      character*(*) c_value,c_param
      integer i_value_chars,ivar
      character*8 form

! Read an integer value (ivar) from the character string (c_value).
      write (form,90) i_value_chars
   90 format ('(i',i2,')')
      read (c_value,form) ivar

!KRA
!      print *,c_param,' =',ivar

      return
      end subroutine char2int

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine get_param_data(input_string,c_param,c_value,&
     &  i_param_chars,i_value_chars,icomment_flag)

      implicit none

      character*(*) input_string,c_param,c_value

      integer leading_blanks,trailing_blanks
      integer i,icomment_flag
      integer i_param_start,i_equals_position,i_param_end,&
     &  i_value_start,i_value_end,i_param_chars,i_value_chars,&
     &  i_loc,i_leading_blanks,i_trailing_blanks

! If the input string is not a comment line, process the data.
      if (input_string(1:1).ne.'!') then

! First count the number of leading and trailing blanks.
        i_leading_blanks = leading_blanks(input_string)
        i_trailing_blanks = trailing_blanks(input_string)

! If the input string is not completely blank, process the data.
        if (i_leading_blanks.ne.len(input_string)) then
          icomment_flag = 0

! Define the starting and ending points of the parameter name and
!   parameter value.
          i_param_start = i_leading_blanks + 1
          i_equals_position = index(input_string,'=')
          i_param_end = i_equals_position - 2
          i_value_start = i_equals_position + 2
          i_value_end =  len(input_string) - i_trailing_blanks
          i_param_chars = i_param_end - i_param_start + 1
          i_value_chars = i_value_end - i_value_start + 1

! Pull out the parameter name and value.
          do i=1,i_param_chars
            i_loc = i + i_param_start - 1
            c_param(i:i) = input_string(i_loc:i_loc)
          enddo

          do i=1,i_value_chars
            i_loc = i + i_value_start - 1
            c_value(i:i) = input_string(i_loc:i_loc)
          enddo

! Check to see whether we have found any paths/fnames that are
!   longer than 80 characters.
          if (i_value_chars.gt.80) then
            print *
            print *, input_string
            print *, 'Number of characters found =',i_value_chars
            print *
            print *, 'Found path/fname longer than 80 characters;'
            print *, 'the limit is 80 characters for the string'
            print *, 'starting in the second character position'
            print *, 'after the "=" sign.'
            print *
            stop
          endif

        else
          icomment_flag = 1
        endif

      else
        icomment_flag = 1
      endif

      return
      end subroutine get_param_data

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      integer function leading_blanks(input_string)

      implicit none

      integer k
      character*(*) input_string

! Count the number of blanks preceeding the first non-blank
!   character.
      leading_blanks = 0
      do k=1,len(input_string)
        if (input_string(k:k).eq.' ') then
          leading_blanks = leading_blanks + 1
        else
          return
        endif
      enddo

      return
      end function leading_blanks

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      integer function trailing_blanks(input_string)

      implicit none

      integer k
      character*(*) input_string

! Count the number of blanks following the last non-blank
!   character.
      trailing_blanks = 0
      do k=len(input_string),1,-1
        if (input_string(k:k).eq.' ') then
          trailing_blanks = trailing_blanks + 1
        else
          return
        endif
      enddo

      return
      end function trailing_blanks


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine GET_DOT_PAR_PATH_LOCATION(snowmodel_dot_par_fname)

! The following provides the ability to read in a 'snowmodel.par'
!   file that is called something other than 'snowmodel.par' and/or
!   located somewhere different than the directory where the
!   'snowmodel' executable file is located.  If you run:
!     snowmodel
!   it assumes the .par file is called 'snowmodel.par' and is
!   located in the default location.  You can also run the code as
!   follows:
!     snowmodel parpath/parname.par
!   and it will look for the .par file in the different location
!   and/or the different name.

      implicit none

      integer nargs,iargc
      character*100 snowmodel_dot_par_fname

      snowmodel_dot_par_fname = "snowmodel.par"

      nargs = iargc()

      if (nargs.eq.1) then
        call getarg (1,snowmodel_dot_par_fname)
      endif

      print *
      print *, snowmodel_dot_par_fname

      return
      end subroutine GET_DOT_PAR_PATH_LOCATION

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

