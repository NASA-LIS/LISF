! preprocess_code.f90

! Perform a variety of preprocessing steps, like read in topography
!   and vegetation arrays, open input and output files, etc.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE PREPROCESS_CODE(topoveg_fname,const_veg_flag,&
     &  vegtype,veg_z0,vegsnowdepth,fetch,xmu,C_z,h_const,&
     &  wind_min,Up_const,dz_susp,ztop_susp,fall_vel,Ur_const,&
     &  ro_water,ro_air,gravity,vonKarman,pi,twopio360,snow_z0,&
     &  nx,ny,sum_sprec,sum_qsubl,sum_trans,sum_unload,topo,&
     &  topo_land,snow_d,topoflag,snow_d_init,snow_d_init_const,&
     &  soft_snow_d,met_input_fname,igrads_metfile,deltax,deltay,&
     &  snowtran_output_fname,micromet_output_fname,&
     &  enbal_output_fname,snowpack_output_fname,print_micromet,&
     &  print_enbal,print_snowpack,print_snowtran,run_micromet,&
     &  run_enbal,run_snowpack,run_snowtran,ro_snow_grid,swe_depth,&
     &  sum_runoff,sum_prec,ro_snow,twolayer_flag,sum_Qcs,&
     &  canopy_int,ascii_topoveg,topo_ascii_fname,icorr_factor_loop,&
     &  veg_ascii_fname,undef,isingle_stn_flag,max_iter,&
     &  i_tair_flag,i_rh_flag,i_wind_flag,i_prec_flag,sum_glacmelt,&
     &  snow_depth,sum_d_canopy_int,corr_factor,icorr_factor_index,&
     &  sum_sfcsublim,barnes_lg_domain,n_stns_used,k_stn,xmn,ymn,&
     &  ro_soft_snow_old,sum_swemelt,xlat,lat_solar_flag,xlat_grid,&
     &  xlon_grid,UTC_flag,dt,swe_depth_old,canopy_int_old,&
     &  vegsnowd_xy,iveg_ht_flag,ihrestart_flag,i_dataassim_loop,&
     &  multilayer_snowpack,max_layers,multilayer_output_fname,&
     &  print_multilayer,KK,tslsnowfall,tsls_threshold,&
     &  irun_data_assim,izero_snow_date,iclear_mn,iclear_dy,&
     &  xclear_hr,snod_layer,swed_layer,ro_layer,T_old,gamma,&
     &  icond_flag,curve_lg_scale_flag,curve_wt_lg,check_met_data,&
     &  seaice_run,snowmodel_line_flag,xg_line,yg_line,print_user,&
     &  cf_precip_flag,cf_precip,print_inc,xhour_init,Tabler_1_flag,&
     &  Tabler_2_flag,iyear_init,imonth_init,iday_init,print_var,&
     &  output_path_wo_assim,output_path_wi_assim,nrecs_max,&
     &  tabler_sfc_path_name,print_outvars,diam_layer)

      use snowmodel_inc
      implicit none

!KRA
      logical snowmodel_masterproc
!KRA

      integer i,j,k,nx,ny,igrads_metfile,n_recs_out,iheader,&
     &  isingle_stn_flag,max_iter,i_tair_flag,i_rh_flag,i_wind_flag,&
     &  i_prec_flag,iter,iobs_num,n_stns_used,nveg,iveg_ht_flag,&
     &  lat_solar_flag,ihrestart_flag,nstns_orig,i_dataassim_loop,&
     &  multilayer_snowpack,max_layers,irun_data_assim,&
     &  izero_snow_date,iclear_mn,iclear_dy,icond_flag,&
     &  iyear_init,imonth_init,iday_init

      real ro_water,ro_air,gravity,vonKarman,snow_z0,&
     &  fetch,xmu,C_z,h_const,wind_min,Up_const,check_met_data,&
     &  dz_susp,ztop_susp,fall_vel,Ur_const,pi,twopio360,topoflag,&
     &  snow_d_init_const,const_veg_flag,ro_snow,twolayer_flag,&
     &  ascii_topoveg,undef,barnes_lg_domain,xlat,UTC_flag,dt,&
     &  print_multilayer,xclear_hr,curve_lg_scale_flag,seaice_run,&
     &  snowmodel_line_flag,print_user,print_inc,xhour_init,&
     &  Tabler_1_flag,Tabler_2_flag

      real topo_land(nx,ny)
      real topo(nx,ny)
      real vegtype(nx,ny)
      real xlat_grid(nx,ny)
      real xlon_grid(nx,ny)

      real snow_d(nx,ny)
      real snow_depth(nx,ny)
      real snow_d_init(nx,ny)
      real canopy_int(nx,ny)
      real swe_depth_old(nx,ny)
      real canopy_int_old(nx,ny)

      real sum_sprec(nx,ny)
      real sum_qsubl(nx,ny)
      real sum_trans(nx,ny)
      real sum_unload(nx,ny)
      real soft_snow_d(nx,ny)
      real ro_soft_snow_old(nx,ny)
      real ro_snow_grid(nx,ny)
      real swe_depth(nx,ny)
      real sum_prec(nx,ny)
      real sum_runoff(nx,ny)
      real sum_Qcs(nx,ny)
      real sum_glacmelt(nx,ny)
      real sum_swemelt(nx,ny)
      real sum_d_canopy_int(nx,ny)
      real sum_sfcsublim(nx,ny)

      real vegsnowdepth(nvegtypes)
      real veg_z0(nx,ny)
      real vegsnowd_xy(nx,ny)

      real curve_wt_lg(nx,ny)

      real corr_factor(nx_max,ny_max,max_obs_dates+1)
      integer icorr_factor_index(max_time_steps)
      integer icorr_factor_loop

      integer k_stn(nx,ny,9)
      double precision xmn,ymn
      double precision nrecs_max,nrecs
      real deltax,deltay
      integer icount,iii,jjj
      double precision xg_line(nx,ny),yg_line(nx,ny)

      real run_micromet,run_enbal,run_snowpack,run_snowtran
      real print_micromet,print_enbal,print_snowpack,print_snowtran

      character*80 topoveg_fname,met_input_fname,topo_ascii_fname,&
     &  veg_ascii_fname
      character*80 snowtran_output_fname,micromet_output_fname,&
     &  enbal_output_fname,snowpack_output_fname,&
     &  multilayer_output_fname

      character*80 tabler_sfc_path_name
      character*80 output_path_wo_assim,output_path_wi_assim

      character*1 print_var(n_print_vars)
      character*4 print_outvars(n_print_vars)

      integer KK(nx,ny)
      real tslsnowfall(nx,ny)
      real tsls_threshold
      real snod_layer(nx,ny,nz_max)
      real swed_layer(nx,ny,nz_max)
      real ro_layer(nx,ny,nz_max)
      real T_old(nx,ny,nz_max)
      real gamma(nx,ny,nz_max)
      real diam_layer(nx,ny,nz_max)

      real cf_precip(nx,ny)
      real cf_precip_flag,cf_precip_scalar

      integer ipath_length,i_len_wo,i_len_wi,trailing_blanks
      character*80 vege_ht_fname

      integer nyears,nyear,nobs_total,nobs_dates,nstns,krec

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! seaice_run = 3.0 is no longer supported.  If seaice_run has been
!   set to 3.0 in the .par file, send a message describing your
!   options.
      if (seaice_run.eq.3.0) then
        print *, 'Eulerian sea ice runs are not longer supported'
        print *, '(seaice_run = 3.0).  If you want to restart this'
        print *, 'option, see the notes in here:'
        print *, '/sm/misc_programs/Eulerian_incremental_remapper/'
        stop
      endif

! Check to see whether the maximum array dimensions defined in
!   snowmodel.inc are adequate for the simulation domain defined
!   in snowmodel.par.
      if (snowmodel_line_flag.eq.0.0) then
        if (seaice_run.eq.3.0) then
          if (nx.ne.nx_max .or. ny.ne.ny_max) then
            print *, 'For a sea ice remapping run, nx==nx_max'
            print *, '  and ny==ny_max in the snowmodel.par and'
            print *, '  snowmodel.inc.'
            stop
          endif
        else
          if (nx+1.gt.nx_max .or. ny+1.gt.ny_max) then
            print *, 'Must increase the value of nx_max or ny_max'
            print *, '  in snowmodel.inc to be greater than nx+1'
            print *, '  and/or ny+1.'
            print *, 'nx_max = ',nx_max,'  ny_max = ',ny_max
            print *, 'nx = ',nx,'  ny = ',ny
            stop
          endif
        endif
      else
        if (nx.ge.nx_max .or. ny.ne.1 .or. ny_max.ne.2) then
          print *, 'For snowmodel_line_flag = 1.0, we suggest setting'
          print *, 'nx = number of grid cells, ny = 1, nx_max = nx+1,'
          print *, 'and ny_max = ny+1 = 2 in snowmodel.inc.'
          print *, '  The current values are:'
          print *, '    nx_max = ',nx_max,'  ny_max = ',ny_max
          print *, '    nx = ',nx,'  ny = ',ny
          stop
        endif
      endif

      if (multilayer_snowpack.eq.1) then
        if (max_layers+1.gt.nz_max) then
          print *, 'nz_max in snowmodel.inc must be at least 1 greater'
          print *, '  than max_layers in the snowmodel.inc file.  So,'
          print *, '  if you want to run the multi-layer snowpack model'
          print *, '  with a single snow layer, set nz_max=2.  If you'
          print *, '  want to run the original single-layer snowpack'
          print *, '  model, you can set nz_max=1 in snowmodel.inc.'
          print *, 'nz_max = ',nz_max
          print *, 'max_layers = ',max_layers
          stop
        endif
      endif

      if (max_iter.gt.max_time_steps) then
        print *, 'Must increase the value of max_time_steps'
        print *, '  in snowmodel.inc to be greater than max_iter.'
        print *, 'max_time_steps = ',max_time_steps
        print *, 'max_iter = ',max_iter
        stop
      endif

! If running the concatenated configuration of the model, check to
!   make sure the rest of the model is configured correctly.
      if (snowmodel_line_flag.eq.1.0) then
        if (run_snowtran.eq.1.0 .and. seaice_run.ne.4.0) then
          print *, 'You cannot run snowmodel_line_flag = 1.0 with'
          print *, 'run_snowtran = 1.0'
        stop
        endif
        if (barnes_lg_domain.eq.0.0) then
          print *, 'If snowmodel_line_flag = 1.0, then you must run'
          print *, 'the model with barnes_lg_domain = 1.0'
        stop
        endif
      endif

! Make sure the time since last snowfall treshold is not less
!   than the model time step.
      if (multilayer_snowpack.eq.1) then
        if (tsls_threshold.lt.dt/3600.0) then
          print *,'Need to correct tsls_threshold to'
          print *,'  be >= dt (in hours).'
          stop
        endif
      endif

! Check the model dt value, and send an error message if dt < 3600.
      if (dt.lt.3600.0) then
        print *, 'You must modify the hour fraction calculation'
        print *, '  in get_model_time subroutine to handle'
        print *, '  dt values less that 3600.0 seconds.'
        print *, 'dt = ',dt
        stop
      endif

! Define the date on which the snow arrays will be zeroed out.
      iclear_mn = izero_snow_date / 10000
      iclear_dy = (izero_snow_date - iclear_mn * 10000) / 100
      xclear_hr = &
     &  real((izero_snow_date - iclear_mn * 10000) - iclear_dy * 100)

! Check to see whether there is enough snow layers to calculate
!   conductive surface fluxes.
      if (icond_flag.eq.1) then
        if (multilayer_snowpack.eq.0 .or. max_layers.lt.2) then
          print *,'If icond_flag = 1, then multilayer_snowpack = 1'
          print *,'  and max_layers >= 2.'
          stop
        endif
      endif

! Read in the topography array.
      if (ascii_topoveg.eq.0.0) then

        open (unit=37,file=topoveg_fname, &
     &    form='unformatted',access='direct',recl=4*nx*ny)
        read (37,rec=1) ((topo_land(i,j),i=1,nx),j=1,ny)

      elseif (ascii_topoveg.eq.1.0) then

! Read off the header lines.  I will assume that all of this
!   information was input in the .par file correctly.
        open (37,file=topo_ascii_fname,form='formatted')
        iheader = 6
        do k=1,iheader
          read (37,*)
        enddo
! Read the data in as real numbers, and do the yrev.
        do j=ny,1,-1
          read (37,*) (topo_land(i,j),i=1,nx)
        enddo

      endif

! If vegetation data is not available on the topography grid,
!   define the vegetation to be constant.
      if (const_veg_flag.ne.0.0) then
        do i=1,nx
          do j=1,ny
            vegtype(i,j) = const_veg_flag
          enddo
        enddo

! Read in the vegetation array.
      else

        if (ascii_topoveg.eq.0.0) then
          read (37,rec=2) ((vegtype(i,j),i=1,nx),j=1,ny)

        elseif (ascii_topoveg.eq.1.0) then

! Read off the header lines.  I will assume that all of this
!   information was input in the .par file correctly.
          open (38,file=veg_ascii_fname,form='formatted')
          do k=1,iheader
            read (38,*)
          enddo
! Read the data in as real numbers, and do the yrev.
          do j=ny,1,-1
            read (38,*) (vegtype(i,j),i=1,nx)
          enddo

        endif

      endif

! Now that we have read in the topo and veg data arrays, check
!   whether all of the values look like valid numbers.
! KRA
      if (ascii_topoveg.eq.0.0 .or. ascii_topoveg.eq.1.0) then
       do i=1,nx
        do j=1,ny
          if (vegtype(i,j).lt.1.0 .or. vegtype(i,j).gt.30.0) then
            print *, 'Found Invalid Vegetation-Type Value'
            print *, '     Value =',vegtype(i,j),'  at ',i,j
            stop
          endif

          if (topo_land(i,j).lt.0.0 .or. topo_land(i,j).gt.9000.0) then
            print *, 'Found Invalid Topography Value'
            print *, '     Value =',topo_land(i,j),'  at ',i,j
            stop
          endif
        enddo
       enddo
! KRA
      endif

! Fill the the vegetation snow-holding depth array for vegetation
!   types 1 through 24 (types 25 through 30 were filled from the
!   .par file.
      call fill_veg_shd(nvegtypes,vegsnowdepth)

! Use vegsnowdepth to fill the 2D spatial array, or read in the
!   user-provided file of the vegetation heights (in m).
      if (iveg_ht_flag.eq.-1) then

! Find the last occurance of '/' in the topo_vege path.
        ipath_length = 0
        do k=1,len(topoveg_fname)
          if (topoveg_fname(k:k).eq.'/') then
            ipath_length = k
          endif
        enddo

        if (ipath_length.eq.0) then
          print *,'vege_ht_fname path not found; looking'
          print *,'  in the SnowModel run directory'
        endif

        vege_ht_fname = &
     &    topoveg_fname(1:ipath_length)//'veg_ht.gdat'

        open (191,file=vege_ht_fname, &
     &    form='unformatted',access='direct',recl=4*nx*ny)
        read (191,rec=1) ((vegsnowd_xy(i,j),i=1,nx),j=1,ny)
        close(191)

      elseif (iveg_ht_flag.eq.1) then

! Find the last occurance of '/' in the veg_ascii_fname path.
        ipath_length = 0
        do k=1,len(veg_ascii_fname)
          if (veg_ascii_fname(k:k).eq.'/') then
            ipath_length = k
          endif
        enddo

        if (ipath_length.eq.0) then
          print *,'vege_ht_fname path not found; looking'
          print *,'  in the SnowModel run directory'
        endif

        vege_ht_fname = &
     &    veg_ascii_fname(1:ipath_length)//'veg_ht.asc'

        open (191,file=vege_ht_fname,form='formatted')
        iheader = 6
        do k=1,iheader
          read (191,*)
        enddo
! Read the data in as real numbers, and do the yrev.
        do j=ny,1,-1
          read (191,*) (vegsnowd_xy(i,j),i=1,nx)
        enddo
        close(191)

      elseif (iveg_ht_flag.eq.0) then

        do i=1,nx
          do j=1,ny
            nveg = nint(vegtype(i,j))
            vegsnowd_xy(i,j) = vegsnowdepth(nveg)
          enddo
        enddo

      endif

! Define the roughness lengths for each of the vegetation types.
!   Note that this is pretty arbitrary, because these values are
!   really only used when there is no blowing snow, and thus have
!   no impact on the simulation except to provide a non-zero value
!   for any such parts of the domain.
      do i=1,nx
        do j=1,ny
          veg_z0(i,j) = 0.25 * vegsnowd_xy(i,j)
        enddo
      enddo

! Read in the large-scale curvature weighting array, if the run
!   requires it.
      if (curve_lg_scale_flag.eq.1.0) then
        open (444,file='extra_met/large_curvature_wt.gdat', &
     &    form='unformatted',access='direct',recl=4*nx*ny)
        read (444,rec=1) ((curve_wt_lg(i,j),i=1,nx),j=1,ny)
        close (444)
      endif

! If this is a sea ice run, open the sea ice concentration file.
      if (seaice_run.ne.0.0) then
        open (445,file='seaice/ice_conc.gdat', &
     &    form='unformatted',access='direct',recl=4*nx*ny)
      endif

! If this is a Lagrangian sea ice parcel trajectory simulation,
!   do some setup checking.
      if (seaice_run.eq.4.0) then
        if (barnes_lg_domain.ne.1.0 .or. n_stns_used.ne.1.0 .or. &
     &    snowmodel_line_flag.ne.1.0) then
          print *
          print *,'if seaice_run = 4.0, then'
          print *,'  barnes_lg_domain = 1.0'
          print *,'  n_stns_used = 1'
          print *,'  snowmodel_line_flag = 1.0'
          print *
          stop
        endif
      endif

! Check to make sure that if you are running SnowTran-3D and the
!   EnBal and SnowPack models, you have also set the flag to run
!   the SnowTran-3D two-layer submodel.
      if (run_enbal.eq.1.0 .and. run_snowpack.eq.1.0 .and. &
     &  run_snowtran.eq.1.0) then
        if (twolayer_flag.ne.1.0) then
          print *, 'For SnowTran-3D with EnBal and SnowPack,'
          print *, '  twolayer_flag must equal 1.0'
          stop
        endif
      endif

! Check to see that the defined domain is large enough to be
!   running SnowTran-3D.
      if (run_snowtran.eq.1.0 .and. seaice_run.ne.4.0) then
        if (nx.lt.3 .or. ny.lt.3) then
          print *, 'To run SnowTran-3D, nx and ny must both be 3'
          print *, '  or greater (see SnowTran-3D code/notes)'
          stop
        endif
      endif

! Check to see whether the model is configured correctly to be
!   running the multi-layer snow model.
      if (multilayer_snowpack.eq.1) then
        if (ihrestart_flag.ne.-2 .or. snow_d_init_const.ne.0.0) then
          print *, 'The multi-layer snowpack model requires:'
          print *, '  ihrestart_flag = -2'
          print *, '  snow_d_init_const = 0.0'
          stop
        endif
      endif

! Get a collection of constants that are not generally modified.
      call constants(fetch,xmu,C_z,h_const,wind_min,Up_const, &
     &  dz_susp,ztop_susp,fall_vel,Ur_const,ro_water,ro_air, &
     &  gravity,vonKarman,pi,twopio360,snow_z0)

! Run a check to see if SnowTran-3D is being run with grid
!   increments that are too large.
      if (deltax.gt.500.0 .or. deltay.gt.500.0) then
        if (seaice_run.eq.0.0) then
          if (run_snowtran.eq.1.0) then
            print *
            print *, '!!! deltax,y should not be greater than 500 m'
            print *, '!!!    if you are also running SnowTran-3D'
            print *
            stop
          endif
        endif
      endif

! Initialize the summing arrays, and define the initial snow-depth
!   distributions.
      call initialize(nx,ny,sum_sprec,sum_qsubl,sum_trans,&
     &  sum_unload,topo,topo_land,snow_d,topoflag,snow_d_init,&
     &  snow_d_init_const,soft_snow_d,ro_water,sum_sfcsublim,&
     &  ro_snow_grid,swe_depth,sum_runoff,sum_prec,ro_snow,&
     &  sum_Qcs,canopy_int,sum_glacmelt,snow_depth,sum_d_canopy_int,&
     &  ro_soft_snow_old,sum_swemelt,swe_depth_old,canopy_int_old,&
     &  ihrestart_flag,i_dataassim_loop,max_iter,corr_factor,&
     &  icorr_factor_index,KK,tslsnowfall,tsls_threshold,snod_layer,&
     &  swed_layer,ro_layer,T_old,gamma,diam_layer)

! Check to see whether the data assimilation has been configured
!   correctly.
      if (irun_data_assim.eq.1) then

! Check to see whether the required output files will be created.
        if (print_user.ne.1.0) then
          print *, 'For a data assimilation run print_user must = 1.0'
          stop
        endif

        if (print_var(19).ne.'y') then
          print *, 'print_var_19 == y for a data assimilation run'
          stop
        endif

        if (print_var(20).ne.'y') then
          print *, 'print_var_20 == y for a data assimilation run'
          stop
        endif

! Check to see whether the corr_factor array is defined in the
!   snowmodel.inc file to be large enough to do the assimilation.
!   max_obs_dates is used in the data assimilation routines.  It
!   must be greater than or equal to the number of observation
!   dates in the entire simulation + (plus) the number of years
!   in the simulation.  For example, for a 6-year simulation with
!   2 observation dates in each year, you would set max_obs_dates
!   to be = 2obs*6yrs+6yrs = 18 or greater.  For a 6-year run with
!   4 observation dates in 2 of the years, and 0 observation dates
!   in the other 4 years, max_obs_dates = 8obs+6yrs = 14 or
!   greater.
!     parameter (max_obs_dates=18)
        if (icorr_factor_loop.eq.1) then
          open (unit=61,file='swe_assim/swe_obs.dat')
          read(61,*) nyears
          nobs_total = 0
          do nyear=1,nyears
            read(61,*) nobs_dates
            if (nobs_dates.gt.0) then
              nobs_total = nobs_total + nobs_dates
              do iobs_num=1,nobs_dates
!               read(61,*) iiyr,iimo,iidy
                read(61,*)
                read(61,*) nstns
                do k=1,nstns
!                 read(61,*) obsid(k),xstn(k),ystn(k),swe_obs(k)
                  read(61,*)
                enddo
              enddo
            endif
          enddo
          close (61)
          krec = nobs_total + nyears
          if (krec.gt.max_obs_dates) then
          print *
          print *, 'For a DATA ASSIMILATION RUN, MAX_OBS_DATES must be'
          print *, 'defined in SNOWMODEL.INC to be greater than the'
          print *, 'number of obs dates in the entire simulation +'
          print *, '(plus) the number of years in the simulation.  For'
          print *, 'example, for a 6-year simulation with 2 observation'
          print *, 'dates in each year, you would set max_obs_dates to'
          print *, 'be = 2obs*6yrs+6yrs = 18 or greater.  For a 6-year'
          print *, 'run with 4 observation dates in 2 of the years,'
          print *, 'and 0 observation dates in the other 4 years,'
          print *, 'max_obs_dates = 8obs+6yrs = 14 or greater.'
          print *
          print *, 'max_obs_dates must be increased in snowmodel.inc'
          print *, 'It looks like you should set max_obs_dates = ',krec
          print *, 'Right now, max_obs_dates = ',max_obs_dates
          print *
          stop
          endif
        endif
      endif

! Initialize the precipitation factor for the first iteration to
!   equal 1.0.
      if (icorr_factor_loop.eq.1) then
        do iobs_num=1,max_obs_dates+1
          do j=1,ny
            do i=1,nx
              corr_factor(i,j,iobs_num) = 1.0
            enddo
          enddo
        enddo
        do iter=1,max_iter
          icorr_factor_index(iter) = 1
        enddo
      endif

! Read or build the latitude array that will be used to do the
!   latitude weighting when calculating incoming solar radiation.
      if (lat_solar_flag.eq.-1) then

        write(*,*) "SM: Reading in gridded latitude binary file"
        open (91,file='extra_met/grid_lat.gdat', &
     &    form='unformatted',access='direct',recl=4*nx*ny)
        read (91,rec=1) ((xlat_grid(i,j),i=1,nx),j=1,ny)
        close(91)

      elseif (lat_solar_flag.eq.1) then

        open (91,file='extra_met/grid_lat.asc',form='formatted')
        iheader = 6
        do k=1,iheader
          read (91,*)
        enddo
! Read the data in as real numbers, and do the yrev.
        do j=ny,1,-1
          read (91,*) (xlat_grid(i,j),i=1,nx)
        enddo
        close(91)

      elseif (lat_solar_flag.eq.0) then

! Print an error if the y-domain is big enough to have important
!   solar radiation differences from south to north.
        if (ny*deltay.gt.500000.0) then
          print *
          print *,'YOUR DOMAIN IS PRETTY BIG TO NOT ACCOUNT FOR'
          print *,'  SOLAR RADIATION VARIATIONS WITH LATITUDE'
          print *,' see the "lat_solar_flag" in snowmodel.par'
          print *
          stop
        endif

        do i=1,nx
          do j=1,ny
            xlat_grid(i,j) = xlat
          enddo
        enddo

      endif

! Read or build the longitude array that will be used to do the
!   longitude influence when calculating incoming solar radiation.
      if (UTC_flag.eq.-1.0) then

        write(*,*) "SM: Reading in gridded longitude binary file"
        open (91,file='extra_met/grid_lon.gdat', &
     &    form='unformatted',access='direct',recl=4*nx*ny)
        read (91,rec=1) ((xlon_grid(i,j),i=1,nx),j=1,ny)
        close(91)

      elseif (UTC_flag.eq.1.0) then

        open (91,file='extra_met/grid_lon.asc',form='formatted')
        iheader = 6
        do k=1,iheader
          read (91,*)
        enddo
! Read the data in as real numbers, and do the yrev.
        do j=ny,1,-1
          read (91,*) (xlon_grid(i,j),i=1,nx)
        enddo
        close(91)

      elseif (UTC_flag.eq.0.0) then

! Print an error if the x-domain is big enough to have important
!   solar radiation differences from east to west.
        if (nx*deltax.gt.500000.0 .and. seaice_run.ne.4) then
          print *
          print *,'YOUR DOMAIN IS PRETTY BIG TO NOT ACCOUNT FOR'
          print *,'  SOLAR RADIATION VARIATIONS WITH LONGITUDE'
          print *,'    see the "UTC_flag" in snowmodel.par'
          print *
        endif

      endif

! Open the MicroMet station data input file.
      if (igrads_metfile.eq.1) then
        open(20,file=met_input_fname,form='unformatted', &
     &    access='direct',recl=4*13)
      else
        open (20,file=met_input_fname,form='formatted')
      endif

! Run a check to see whether there are any time slices with no
!   valid data.
      if (check_met_data.eq.1.0) then
        print *
        print *,'Checking for sufficient met forcing data to'
        print *,'  complete the model simulation.  This may'
        print *,'  take a while, depending on how big your met'
        print *,'  input file is.'
        print *
        call met_data_check(undef,isingle_stn_flag,igrads_metfile, &
     &    max_iter,i_tair_flag,i_rh_flag,i_wind_flag,i_prec_flag)
      endif

! If the concatenated configuration of the model is used, read
!   in the x and y coordinates for the concatenated grid cells.
      if (snowmodel_line_flag.eq.1.0) then
        open (1331,file='extra_met/snowmodel_line_pts.dat')
        do j=1,ny
          do i=1,nx
            read (1331,*) icount,iii,jjj,xg_line(i,j),yg_line(i,j)
          enddo
        enddo
        close (1331)
      endif

! If the large-domain barnes oi scheme is used, generate the
!   nearest-station indexing array.
      if (barnes_lg_domain.eq.1.0) then
!      if (barnes_lg_domain.eq.1.0 .and. snowmodel_masterproc) then
        print *
        print *,'You are running the large-domain Barnes oi scheme'
        print *,'  This requires:'
        print *,'  1) no missing data for the fields of interest'
        print *,'  2) no missing stations during the simulation' 
        print *,'  3) met file must list stations in the same order'
        print *,'  4) the number of nearest stations used is 9 or less'
        print *,'  5)  **** no error checking for this is done ****'
        print *
        print *,'Generating nearest-station index.  Be patient.'
        print *
        if (n_stns_used.gt.9 .or. n_stns_used.lt.1) then
          print *,'invalid n_stns_used value'
          stop
        endif
! KRA
!        call get_nearest_stns_1(nx,ny,xmn,ymn,deltax,deltay,
!     &    n_stns_used,k_stn,snowmodel_line_flag,xg_line,yg_line)
! KRA
      endif

! If this is a history restart run, advance the micromet input
!   file to the restart time.
      if (ihrestart_flag.ge.0) then
        if (igrads_metfile.eq.0) then
          do iter=1,ihrestart_flag
            if (isingle_stn_flag.eq.1) then
              nstns_orig = 1
            else
              read(20,*) nstns_orig
            endif
            do k=1,nstns_orig
              read(20,*)
            enddo
          enddo
        endif
      endif

! Open the files to be used to store model output.

! nrecs_max corresponds to the approximately 2.1 GB Fortran
!   array limit for direct access binary inputs and outputs.
!   The number listed here corresponds to nx = ny = 23170.
!   "nrecs_max * 4 bytes per number" gives the 2.1 GB limit.
      nrecs_max = 536848900

!   For MicroMet.
      if (run_micromet.eq.1.0 .and. print_micromet.eq.1.0) then
        n_recs_out = 9
        nrecs = n_recs_out * nx * ny
        if (nrecs.gt.nrecs_max) then
          print *,'Your simulation domain has too many grid cells'
          print *,'to print the micromet.gdat file.  You must set'
          print *,'print_micromet = 0.0 and use print_user = 1.0.'
          stop
        else
          if (icorr_factor_loop.eq.2) close (81)
          open (81,file=micromet_output_fname,&
     &      form='unformatted',access='direct',recl=4*n_recs_out*nx*ny,&
     &      status='replace')
        endif
      endif

!   For EnBal.
      if (run_enbal.eq.1.0 .and. print_enbal.eq.1.0) then
        n_recs_out = 11
        nrecs = n_recs_out * nx * ny
        if (nrecs.gt.nrecs_max) then
          print *,'Your simulation domain has too many grid cells'
          print *,'to print the enbal.gdat file.  You must set'
          print *,'print_enbal = 0.0 and use print_user = 1.0.'
          stop
        else
          if (icorr_factor_loop.eq.2) close (82)
          open (82,file=enbal_output_fname,&
     &      form='unformatted',access='direct',recl=4*n_recs_out*nx*ny,&
     &      status='replace')
        endif
      endif

!   For SnowPack.
      if (run_snowpack.eq.1.0 .and. print_snowpack.eq.1.0) then
        n_recs_out = 16
        nrecs = n_recs_out * nx * ny
        if (nrecs.gt.nrecs_max) then
          print *,'Your simulation domain has too many grid cells'
          print *,'to print the snowpack.gdat file.  You must set'
          print *,'print_snowpack = 0.0 and use print_user = 1.0.'
          stop
        else
          if (icorr_factor_loop.eq.2) close (83)
          open (83,file=snowpack_output_fname, &
     &      form='unformatted',access='direct',recl=4*n_recs_out*nx*ny,&
     &      status='replace')
        endif
      endif

!   For SnowTran-3D.
      if (run_snowtran.eq.1.0 .and. print_snowtran.eq.1.0) then
        n_recs_out = 7
        nrecs = n_recs_out * nx * ny
        if (nrecs.gt.nrecs_max) then
          print *,'Your simulation domain has too many grid cells'
          print *,'to print the snowtran.gdat file.  You must set'
          print *,'print_snowtran = 0.0 and use print_user = 1.0.'
          stop
        else
          if (icorr_factor_loop.eq.2) close (84)
          open (84,file=snowtran_output_fname,&
     &      form='unformatted',access='direct',recl=4*n_recs_out*nx*ny,&
     &      status='replace')
        endif
      endif

!   For Multi-Layer SnowPack.
      if (run_snowpack.eq.1.0 .and. multilayer_snowpack.eq.1 .and.&
     &  print_multilayer.eq.1.0) then
        nrecs = 4 * nx * ny + 4 * nx * ny * nz_max
        if (nrecs.gt.nrecs_max) then
          print *,'Your simulation domain has too many grid cells'
          print *,'to print the multilayer.gdat file.  Since you'
          print *,'clearly want this information, it will be written'
          print *,'to the directory you defined in the .par file for'
          print *,'the parameter "output_path_wo_assim".'
        else
          if (icorr_factor_loop.eq.2) close (401)
          open (401,file=multilayer_output_fname, &
     &      form='unformatted',access='direct', &
     &      recl=4*(4*nx*ny+4*nx*ny*nz_max), &
     &      status='replace')
        endif
      elseif (run_snowpack.eq.1.0 .and. multilayer_snowpack.eq.1 .and. &
     &  print_multilayer.eq.2.0) then
        nrecs = nx * ny * nz_max
        if (nrecs.gt.nrecs_max) then
          print *,'Your simulation domain has too many grid cells'
          print *,'to print the print_multilayer = 2.0 files.  You'
          print *,'will have to restructure the write statements.'
          print *,'See the example in the outputs_user.f subroutine'
          print *,'where it does the "if (nrecs.gt.nrecs_max) then"'
          print *,'test.'
          stop
        else

         if (icorr_factor_loop.eq.1) then

          i_len_wo = 80 - trailing_blanks(output_path_wo_assim)
          open (401,&
     &    file=output_path_wo_assim(1:i_len_wo)//'multilayer_2Dxy.gdat',&
     &      form='unformatted',access='direct',&
     &      recl=4*(4*nx*ny),status='replace')
          open (402,&
     &    file=output_path_wo_assim(1:i_len_wo)//'multilayer_snod.gdat',&
     &      form='unformatted',access='direct',&
     &      recl=4*nx*ny*nz_max,status='replace')
          open (403,&
     &    file=output_path_wo_assim(1:i_len_wo)//'multilayer_sden.gdat',&
     &      form='unformatted',access='direct',&
     &      recl=4*nx*ny*nz_max,status='replace')
          open (404,&
     &    file=output_path_wo_assim(1:i_len_wo)//'multilayer_swed.gdat',&
     &      form='unformatted',access='direct',&
     &      recl=4*nx*ny*nz_max,status='replace')
          open (405,&
     &    file=output_path_wo_assim(1:i_len_wo)//'multilayer_diam.gdat',&
     &      form='unformatted',access='direct',&
     &      recl=4*nx*ny*nz_max,status='replace')
          open (406,&
     &    file=output_path_wo_assim(1:i_len_wo)//'multilayer_flux.gdat',&
     &      form='unformatted',access='direct',&
     &      recl=4*nx*ny*nz_max,status='replace')
          open (407,&
     &    file=output_path_wo_assim(1:i_len_wo)//'multilayer_temp.gdat',&
     &      form='unformatted',access='direct',&
     &      recl=4*nx*ny*nz_max,status='replace')
          open (408,&
     &    file=output_path_wo_assim(1:i_len_wo)//'multilayer_cond.gdat',&
     &      form='unformatted',access='direct',&
     &      recl=4*nx*ny*nz_max,status='replace')

         elseif (icorr_factor_loop.eq.2) then

          close (401)
          close (402)
          close (403)
          close (404)
          close (405)
          close (406)
          close (407)
          close (408)

          i_len_wi = 80 - trailing_blanks(output_path_wi_assim)
          open (401,&
     &    file=output_path_wi_assim(1:i_len_wi)//'multilayer_2Dxy.gdat',&
     &      form='unformatted',access='direct',&
     &      recl=4*(4*nx*ny),status='replace')
          open (402,&
     &    file=output_path_wi_assim(1:i_len_wi)//'multilayer_snod.gdat',&
     &      form='unformatted',access='direct',&
     &      recl=4*nx*ny*nz_max,status='replace')
          open (403,&
     &    file=output_path_wi_assim(1:i_len_wi)//'multilayer_sden.gdat',&
     &      form='unformatted',access='direct',&
     &      recl=4*nx*ny*nz_max,status='replace')
          open (404,&
     &    file=output_path_wi_assim(1:i_len_wi)//'multilayer_swed.gdat',&
     &      form='unformatted',access='direct',&
     &      recl=4*nx*ny*nz_max,status='replace')
          open (405,&
     &    file=output_path_wi_assim(1:i_len_wi)//'multilayer_diam.gdat',&
     &      form='unformatted',access='direct',&
     &      recl=4*nx*ny*nz_max,status='replace')
          open (406,&
     &    file=output_path_wi_assim(1:i_len_wi)//'multilayer_flux.gdat',&
     &      form='unformatted',access='direct',&
     &      recl=4*nx*ny*nz_max,status='replace')
          open (407,&
     &    file=output_path_wi_assim(1:i_len_wi)//'multilayer_temp.gdat',&
     &      form='unformatted',access='direct',&
     &      recl=4*nx*ny*nz_max,status='replace')
          open (408,&
     &    file=output_path_wi_assim(1:i_len_wi)//'multilayer_cond.gdat',&
     &      form='unformatted',access='direct',&
     &      recl=4*nx*ny*nz_max,status='replace')

         endif
        endif
      endif

! Read in the precipitation correction factor array.
      if (cf_precip_flag.eq.1.0) then

        open (unit=144,file='precip_cf/cf_precip.gdat',&
     &    form='unformatted',access='direct',recl=4*nx*ny)
        read (144,rec=1) ((cf_precip(i,j),i=1,nx),j=1,ny)

      elseif (cf_precip_flag.eq.2.0) then

! Read off the header lines.  I will assume that all of this
!   information was input in the .par file correctly.
        open (144,file='precip_cf/cf_precip.asc',form='formatted')
        iheader = 6
        do k=1,iheader
          read (144,*)
        enddo
! Read the data in as real numbers, and do the yrev.
        do j=ny,1,-1
          read (144,*) (cf_precip(i,j),i=1,nx)
        enddo

      elseif (cf_precip_flag.eq.3.0) then

        open (144,file='precip_cf/cf_precip.dat',form='formatted')
        read (144,*) cf_precip_scalar
        do j=1,ny
          do i=1,nx
            cf_precip(i,j) = cf_precip_scalar
          enddo
        enddo

      endif

! This must be closed so it can be reread if there is a (second)
!   data assimilation loop.
      close (144)

! Generate all of the GrADS control (.ctl) files that correspond
!   to all of the GrADS output (.gdat) files that were generated as
!   part of this model run.  They are (mostly) all placed in a
!   directory called "ctl_files".
      call mk_ctl_files(nx,ny,deltax,deltay,xmn,ymn,dt,&
     &  print_inc,iyear_init,imonth_init,iday_init,xhour_init,&
     &  max_iter,undef,output_path_wo_assim,output_path_wi_assim,&
     &  print_micromet,micromet_output_fname,print_enbal,&
     &  enbal_output_fname,print_snowpack,snowpack_output_fname,&
     &  print_snowtran,snowtran_output_fname,Tabler_1_flag,&
     &  tabler_sfc_path_name,Tabler_2_flag,irun_data_assim,&
     &  print_var,print_outvars,print_multilayer,&
     &  multilayer_output_fname)

! If this is going to be a SnowTran-3D run, print the Copyright
!   header.
      if (run_snowtran.eq.1.0) then
        print *
        print *,&
     & 'cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc'
        print *,&
     & 'c     Snow-Transport Modeling System - 3D (SnowTran-3D)    c'
        print *,&
     & 'c                    Copyright (C) 1998                    c'
        print *,&
     & 'c          by Glen E. Liston, InterWorks Consulting        c'
        print *,&
     & 'c                    All Rights Reserved                   c'
        print *,&
     & 'cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc'
        print *
        print *
      endif

      return
      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine initialize(nx,ny,sum_sprec,sum_qsubl,sum_trans,&
     &  sum_unload,topo,topo_land,snow_d,topoflag,snow_d_init,&
     &  snow_d_init_const,soft_snow_d,ro_water,sum_sfcsublim,&
     &  ro_snow_grid,swe_depth,sum_runoff,sum_prec,ro_snow,&
     &  sum_Qcs,canopy_int,sum_glacmelt,snow_depth,sum_d_canopy_int,&
     &  ro_soft_snow_old,sum_swemelt,swe_depth_old,canopy_int_old,&
     &  ihrestart_flag,i_dataassim_loop,max_iter,corr_factor,&
     &  icorr_factor_index,KK,tslsnowfall,tsls_threshold,snod_layer,&
     &  swed_layer,ro_layer,T_old,gamma,diam_layer)

      use snowmodel_inc
      implicit none

      integer i,j,nx,ny,ihrestart_flag,i_dataassim_loop,max_iter,k

      real topoflag,snow_d_init_const,ro_water,ro_snow
      real sum_sprec(nx,ny)
      real sum_qsubl(nx,ny)
      real sum_trans(nx,ny)
      real sum_unload(nx,ny)
      real topo(nx,ny)
      real topo_land(nx,ny)
      real snow_d(nx,ny)
      real snow_depth(nx,ny)
      real soft_snow_d(nx,ny)
      real ro_soft_snow_old(nx,ny)
      real ro_snow_grid(nx,ny)
      real swe_depth(nx,ny)
      real sum_prec(nx,ny)
      real sum_runoff(nx,ny)
      real sum_Qcs(nx,ny)
      real canopy_int(nx,ny)
      real sum_glacmelt(nx,ny)
      real sum_swemelt(nx,ny)
      real sum_d_canopy_int(nx,ny)
      real sum_sfcsublim(nx,ny)
      real snow_d_init(nx,ny)
      real swe_depth_old(nx,ny)
      real canopy_int_old(nx,ny)

      integer KK(nx,ny)
      real tslsnowfall(nx,ny)
      real tsls_threshold
      real snod_layer(nx,ny,nz_max)
      real swed_layer(nx,ny,nz_max)
      real ro_layer(nx,ny,nz_max)
      real T_old(nx,ny,nz_max)
      real gamma(nx,ny,nz_max)
      real diam_layer(nx,ny,nz_max)

      integer icorr_factor_index(max_time_steps)
      real corr_factor(nx_max,ny_max,max_obs_dates+1)

      if (ihrestart_flag.ge.0) then

! Read in the saved data.
        CALL HRESTART_READ(nx,ny,snow_d,snow_depth,&
     &    canopy_int,soft_snow_d,ro_snow_grid,swe_depth,&
     &    ro_soft_snow_old,snow_d_init,swe_depth_old,&
     &    canopy_int_old,topo,sum_sprec,ihrestart_flag,&
     &    i_dataassim_loop)

        if (i_dataassim_loop.lt.0.0) then
          CALL HRESTART_READ_DA(nx,ny,max_iter,corr_factor,&
     &      icorr_factor_index,i_dataassim_loop)
        endif

        do i=1,nx
          do j=1,ny
! Fill the summing arrays.
            sum_runoff(i,j) = 0.0
            sum_prec(i,j) = 0.0
!           sum_sprec(i,j) = 0.0
            sum_qsubl(i,j) = 0.0
            sum_trans(i,j) = 0.0
            sum_unload(i,j) = 0.0
            sum_Qcs(i,j) = 0.0
            sum_glacmelt(i,j) = 0.0
            sum_swemelt(i,j) = 0.0
            sum_d_canopy_int(i,j) = 0.0
            sum_sfcsublim(i,j) = 0.0

! Define the initial snow-depth distributions.
!           snow_d_init(i,j) = snow_d_init_const
!           snow_d(i,j) = snow_d_init(i,j)
!           snow_depth(i,j) = snow_d_init(i,j)
!           canopy_int(i,j) = 0.0
!           soft_snow_d(i,j) = snow_d(i,j)
!           ro_snow_grid(i,j) = ro_snow
!           swe_depth(i,j) = snow_d(i,j) * ro_snow_grid(i,j) / ro_water
!           ro_soft_snow_old(i,j) = 50.0
!           swe_depth_old(i,j) = swe_depth(i,j)
!           canopy_int_old(i,j) = canopy_int(i,j)
          enddo
        enddo

      else

        do i=1,nx
          do j=1,ny
! Fill the summing arrays.
            sum_runoff(i,j) = 0.0
            sum_prec(i,j) = 0.0
            sum_sprec(i,j) = 0.0
            sum_qsubl(i,j) = 0.0
            sum_trans(i,j) = 0.0
            sum_unload(i,j) = 0.0
            sum_Qcs(i,j) = 0.0
            sum_glacmelt(i,j) = 0.0
            sum_swemelt(i,j) = 0.0
            sum_d_canopy_int(i,j) = 0.0
            sum_sfcsublim(i,j) = 0.0

! Define the initial snow-depth distributions.
            snow_d_init(i,j) = snow_d_init_const
            snow_d(i,j) = snow_d_init(i,j)
            snow_depth(i,j) = snow_d_init(i,j)
            canopy_int(i,j) = 0.0
            soft_snow_d(i,j) = snow_d(i,j)
            ro_snow_grid(i,j) = ro_snow
            swe_depth(i,j) = snow_d(i,j) * ro_snow_grid(i,j) / ro_water
            ro_soft_snow_old(i,j) = 50.0
            swe_depth_old(i,j) = swe_depth(i,j)
            canopy_int_old(i,j) = canopy_int(i,j)

! Initialize the multi-layer snowpack arrays.
            KK(i,j) = 0
            tslsnowfall(i,j) = tsls_threshold
          enddo
        enddo

        do i=1,nx
          do j=1,ny
            do k=1,nz_max
              snod_layer(i,j,k) = 0.0
              swed_layer(i,j,k) = 0.0
              ro_layer(i,j,k) = ro_snow
              T_old(i,j,k) = 273.15
              gamma(i,j,k) = 0.138 - 1.01 * (ro_layer(i,j,k)/1000.0) + &
     &          3.233 * (ro_layer(i,j,k)/1000.0)**2
              diam_layer(i,j,k) = 0.5 / 1000.0
            enddo
          enddo
        enddo

        if (topoflag.eq.1.0) then
          do i=1,nx
            do j=1,ny
              topo(i,j) = topo_land(i,j) + snow_d(i,j)
            enddo
          enddo
        elseif (topoflag.eq.0.0) then
          do i=1,nx
            do j=1,ny
              topo(i,j) = topo_land(i,j)
            enddo
          enddo
        endif

      endif

      return
      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine constants(fetch,xmu,C_z,h_const,wind_min,Up_const,&
     &  dz_susp,ztop_susp,fall_vel,Ur_const,ro_water,ro_air,&
     &  gravity,vonKarman,pi,twopio360,snow_z0)

      implicit none

      real fetch,xmu,C_z,h_const,wind_min,Up_const,&
     &  dz_susp,ztop_susp,fall_vel,Ur_const,ro_water,ro_air,&
     &  gravity,vonKarman,pi,twopio360,snow_z0

! These constants are not generally modified for a particular model
!   run.

! Snow surface roughness length.
      snow_z0 = 0.001

! Constants related to surface shear stress and saltation
!   transport.
      fetch = 500.0
      xmu = 3.0
      C_z = 0.12
      h_const = 1.6
      wind_min = 4.0

! Constants related to suspended snow profile.
      Up_const = 2.8
      dz_susp = 0.20
      ztop_susp = 2.0
      fall_vel = 0.3
      Ur_const = 0.5

! General constants.
      ro_water = 1000.0
      ro_air = 1.275
      gravity = 9.81
      vonKarman = 0.4
      pi = 2.0 * acos(0.0)
      twopio360 = 2.0 * pi / 360.0

      return
      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine fill_veg_shd(nvegtypes,vegsnowdepth)

      implicit none

      integer k,nvegtypes,nvegtypes_fixed

      parameter (nvegtypes_fixed=24)

      real vegsnowdepth(nvegtypes),vegsnowdepth_fixed(nvegtypes_fixed)

! Fill the the vegetation snow-holding depth array for
!   vegetation types 1 through 24 (types 25 through 30 were filled
!   from the .par file.
!
! The following summary was taken from the .par file.
!
! The vegetation types are assumed to range from 1 through 30.  The
!   last 6 types are available to be user-defined vegetation types
!   and vegetation snow-holding depth.  The first 24 vegetation
!   types, and the associated vegetation snow-holding depth
!   (meters), are hard-coded to be:
!
! code description           veg_shd  example                    class
!
!  1  coniferous forest       15.00  spruce-fir/taiga/lodgepole  forest
!  2  deciduous forest        12.00  aspen forest                forest
!  3  mixed forest            14.00  aspen/spruce-fir/low taiga  forest
!  4  scattered short-conifer  8.00  pinyon-juniper              forest
!  5  clearcut conifer         4.00  stumps and regenerating     forest
! 
!  6  mesic upland shrub       0.50  deeper soils, less rocky    shrub
!  7  xeric upland shrub       0.25  rocky, windblown soils      shrub
!  8  playa shrubland          1.00  greasewood, saltbush        shrub
!  9  shrub wetland/riparian   1.75  willow along streams        shrub
! 10  erect shrub tundra       0.65  arctic shrubland            shrub
! 11  low shrub tundra         0.30  low to medium arctic shrubs shrub
! 
! 12  grassland rangeland      0.15  graminoids and forbs        grass
! 13  subalpine meadow         0.25  meadows below treeline      grass
! 14  tundra (non-tussock)     0.15  alpine, high arctic         grass
! 15  tundra (tussock)         0.20  graminoid and dwarf shrubs  grass
! 16  prostrate shrub tundra   0.10  graminoid dominated         grass
! 17  arctic gram. wetland     0.20  grassy wetlands, wet tundra grass
! 
! 18  bare                     0.01                              bare
!
! 19  water/possibly frozen    0.01                              water
! 20  permanent snow/glacier   0.01                              water
! 
! 21  residential/urban        0.01                              human
! 22  tall crops               0.40  e.g., corn stubble          human
! 23  short crops              0.25  e.g., wheat stubble         human
! 24  ocean                    0.01                              water
!
! 25  user defined (see below)
! 26  user defined (see below)
! 27  user defined (see below)
! 28  user defined (see below)
! 29  user defined (see below)
! 30  user defined (see below)
!
! Define the vegetation snow-holding depth (meters) for each
!   of the user-defined vegetation types.  The numbers in the
!   list below correspond to the vegetation-type numbers in the
!   vegetation-type data array (veg type 25.0 -> veg_shd_25).  Note
!   that even if these are not used, they cannot be commented out
!   or you will get an error message.
!     veg_shd_25 = 0.10
!     veg_shd_26 = 0.10
!     veg_shd_27 = 0.10
!     veg_shd_28 = 0.10
!     veg_shd_29 = 0.10
!     veg_shd_30 = 0.10

      data vegsnowdepth_fixed/15.00, 12.00, 14.00,  8.00,  4.00, &
     &                         0.50,  0.25,  1.00,  1.75,  0.65,  0.30,&
     &                         0.15,  0.25,  0.15,  0.20,  0.10,  0.20,&
     &                         0.01,  0.01,  0.01,  0.01,  0.40,  0.25,&
     &                         0.01/

      do k=1,nvegtypes_fixed
        vegsnowdepth(k) = vegsnowdepth_fixed(k)
      enddo

      return
      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine met_data_check(undef,isingle_stn_flag,igrads_metfile,&
     &  max_iter,i_tair_flag,i_rh_flag,i_wind_flag,i_prec_flag)

      use snowmodel_inc
      implicit none

      integer iyr,imo,idy      ! year, month, and day of data
      real xhr                 ! decimal hour
      integer idstn            ! station id number

      integer k,nstns_orig,isingle_stn_flag,igrads_metfile,iter,&
     &  n_good_Tair,n_good_rh,n_good_wspd,n_good_wdir,n_good_prec,&
     &  n_notgood_vars,i_tair_flag,i_rh_flag,i_wind_flag,i_prec_flag,&
     &  max_iter

      real Tair_orig(nstns_max),rh_orig(nstns_max)
      real winddir_orig(nstns_max),windspd_orig(nstns_max)
      double precision xstn_orig(nstns_max),ystn_orig(nstns_max)
      real elev_orig(nstns_max),prec_orig(nstns_max)
      real undef               ! undefined value
      real elevation_flag

      n_notgood_vars = 0
      elevation_flag = 0.0

      do iter=1,max_iter

        n_good_Tair = 0
        n_good_rh = 0
        n_good_wspd = 0
        n_good_wdir = 0
        n_good_prec = 0

        if (igrads_metfile.eq.1) then
          nstns_orig = 1
        else
          if (isingle_stn_flag.eq.1) then
            nstns_orig = 1
          else
            read(20,*) nstns_orig
          endif
        endif

        if (nstns_orig.gt.nstns_max) then
          print *, 'The number of met stations in your MicroMet'
          print *, 'input file exceeds nstns_max in snowmodel.inc.'
          print *, 'This occurs at iter =',iter
          stop
        endif

        do k=1,nstns_orig

          if (igrads_metfile.eq.1) then
            read(20,rec=iter) iyr,imo,idy,xhr,idstn,xstn_orig(k),&
     &        ystn_orig(k),elev_orig(k),Tair_orig(k),rh_orig(k),&
     &        windspd_orig(k),winddir_orig(k),prec_orig(k)
          else
            read(20,*) iyr,imo,idy,xhr,idstn,xstn_orig(k),&
     &        ystn_orig(k),elev_orig(k),Tair_orig(k),rh_orig(k),&
     &        windspd_orig(k),winddir_orig(k),prec_orig(k)
          endif

! Check for any NaN values.  They are not allowed.
          if (Tair_orig(k).ne.Tair_orig(k) .or.&
     &      rh_orig(k).ne.rh_orig(k) .or.&
     &      windspd_orig(k).ne.windspd_orig(k) .or.&
     &      winddir_orig(k).ne.winddir_orig(k) .or.&
     &      prec_orig(k).ne.prec_orig(k)) then
            print *
            print *,'  YOU HAVE NaN VALUES IN YOUR MET FORCING INPUT'
            print *,'  FILE.  THEY ARE NOT ALLOWED ANYWHERE IN THE'
            print *,'  MicroMet INPUT FILE.  THIS MUST BE CORRECTED'
            print *,'  BEFORE YOU CAN CONTINUE.'
            print *
            stop
          endif

! Count the good values at this time.
          if (Tair_orig(k).ne.undef) n_good_Tair = n_good_Tair + 1
          if (rh_orig(k).ne.undef) n_good_rh = n_good_rh + 1
          if (windspd_orig(k).ne.undef) n_good_wspd = n_good_wspd + 1
          if (winddir_orig(k).ne.undef) n_good_wdir = n_good_wdir + 1
          if (prec_orig(k).ne.undef) n_good_prec = n_good_prec + 1

          if (elev_orig(k).lt.0.0) then
            elevation_flag = 1.0
            print *,'elevation = ',elev_orig(k),'  for stn id = ',idstn
          endif

        enddo

! Check to see whether there are any variables with no valid data
!   at this time slice.
        if (n_good_Tair.eq.0 .and. i_tair_flag.eq.1) then
          n_notgood_vars = n_notgood_vars + 1
          print *,'no good Tair data at           ',iyr,imo,idy,xhr
        endif

        if (n_good_rh.eq.0 .and. i_rh_flag.eq.1) then
          n_notgood_vars = n_notgood_vars + 1
          print *,'no good rh data at             ',iyr,imo,idy,xhr
        endif

        if (n_good_wspd.eq.0 .and. i_wind_flag.eq.1) then
          n_notgood_vars = n_notgood_vars + 1
          print *,'no good wind speed data at     ',iyr,imo,idy,xhr
        endif

        if (n_good_wdir.eq.0 .and. i_wind_flag.eq.1) then
          n_notgood_vars = n_notgood_vars + 1
          print *,'no good wind direction data at ',iyr,imo,idy,xhr
        endif

        if (n_good_prec.eq.0 .and. i_prec_flag.eq.1) then
          n_notgood_vars = n_notgood_vars + 1
          print *,'no good precipitation data at  ',iyr,imo,idy,xhr
        endif

      enddo

      if (n_notgood_vars.gt.0) then
        print *
        print *,' FOUND TIMES WITH NO VALID MET OBSERVATIONS'
        print *,'NEED TO CORRECT THE PROBLEM BEFORE CONTINUING'
        stop
      endif

      if (elevation_flag.eq.1.0) then
        print *
        print *,' FOUND A NEGATIVE OR UNDEFINED STATION ELEVATION.'
        print *,'STATION ELEVATIONS CANNOT BE UNDEFINED, BUT THEY.'
        print *,'CAN BE NEGATIVE FOR A PLACE LIKE DEATH VALLEY.'
        print *,'YOU NEED TO CORRECT ANY PROBLEMS BEFORE CONTINUING.'
        stop
      endif

      if (igrads_metfile.eq.0) rewind (20)

      return
      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine get_nearest_stns_1(nx,ny,xmn,ymn,deltax,deltay,&
     &  n_stns_used,k_stn,snowmodel_line_flag,xg_line,yg_line)

      use snowmodel_inc
      implicit none

      double precision xstn(nstns_max)
      double precision ystn(nstns_max)
      double precision dsq(nstns_max)
      double precision xg_line(nx,ny),yg_line(nx,ny)
      real snowmodel_line_flag

      double precision xg,yg,xmn,ymn,dist_min
      real deltax,deltay,x1,x2,x3,x4,x5,x6,x7

      integer i,j,k,kk,nstns,n_stns_used,nx,ny,i1,i2,i3,i4
      integer k_stn(nx,ny,9)

! Read the station information for the first (and all) time step(s).
      read(20,*) nstns
      do k=1,nstns
        read(20,*) i1,i2,i3,x1,i4,xstn(k),ystn(k),&
     &    x2,x3,x4,x5,x6,x7
      enddo
      rewind (20)

      do j=1,ny
        do i=1,nx

! xcoords of grid nodes at index i,j
! ycoords of grid nodes at index i,j
          if (snowmodel_line_flag.eq.1.0) then
            xg = xg_line(i,j)
            yg = yg_line(i,j)
          else
            xg = xmn + deltax * (real(i) - 1.0)
            yg = ymn + deltay * (real(j) - 1.0)
          endif

! Loop through all of the stations, calculating the distance
!   between the current grid point and each of the stations.
          do k=1,nstns
            dsq(k) = (xg - xstn(k))**2 + (yg - ystn(k))**2
          enddo

! Loop through each of the station distances and find the
!   stations closest to the grid point in question.
          do kk=1,n_stns_used
            dist_min = 1.0e30
            do k=1,nstns
              if (dsq(k).le.dist_min) then
                k_stn(i,j,kk) = k
                dist_min = dsq(k)
              endif
            enddo

! Eliminate the last found minimum from the next search by making
!   its distance a big number.
            dsq(k_stn(i,j,kk)) = 1.0e30
          enddo

        enddo
      enddo

      return
      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine HRESTART_SAVE(nx,ny,iter,snow_d,snow_depth,&
     &  canopy_int,soft_snow_d,ro_snow_grid,swe_depth,&
     &  ro_soft_snow_old,snow_d_init,swe_depth_old,&
     &  canopy_int_old,topo,sum_sprec,icorr_factor_loop,&
     &  max_iter)

      use snowmodel_inc
      implicit none

      integer i,j,nx,ny,iter,icorr_factor_loop,max_iter
      real snow_d(nx,ny)
      real snow_depth(nx,ny)
      real canopy_int(nx,ny)
      real soft_snow_d(nx,ny)
      real ro_snow_grid(nx,ny)
      real swe_depth(nx,ny)
      real ro_soft_snow_old(nx,ny)
      real snow_d_init(nx,ny)
      real swe_depth_old(nx,ny)
      real canopy_int_old(nx,ny)
      real topo(nx,ny)
      real sum_sprec(nx,ny)

      character*18 name1
      character*1 name2
      character*5 name3
      character*5 niter
      character*1 iloop
      character*30 fname

! Build the file name so it includes the interation number.
      name1 = 'hrestart/hrestart_'
      name2 = '_'
      name3 = '.gdat'

      write(niter,'(i5.5)') iter
      write(iloop,'(i1.1)') icorr_factor_loop
      fname = name1//niter//name2//iloop//name3

! Save the data.
      open(151,file=fname,&
     &  form='unformatted',access='direct',recl=4*nx*ny)

      write (151,rec=1) ((snow_d(i,j),i=1,nx),j=1,ny)
      write (151,rec=2) ((snow_depth(i,j),i=1,nx),j=1,ny)
      write (151,rec=3) ((canopy_int(i,j),i=1,nx),j=1,ny)
      write (151,rec=4) ((soft_snow_d(i,j),i=1,nx),j=1,ny)
      write (151,rec=5) ((ro_snow_grid(i,j),i=1,nx),j=1,ny)
      write (151,rec=6) ((swe_depth(i,j),i=1,nx),j=1,ny)
      write (151,rec=7) ((ro_soft_snow_old(i,j),i=1,nx),j=1,ny)
      write (151,rec=8) ((snow_d_init(i,j),i=1,nx),j=1,ny)
      write (151,rec=9) ((swe_depth_old(i,j),i=1,nx),j=1,ny)
      write (151,rec=10) ((canopy_int_old(i,j),i=1,nx),j=1,ny)
      write (151,rec=11) ((topo(i,j),i=1,nx),j=1,ny)
      write (151,rec=12) ((sum_sprec(i,j),i=1,nx),j=1,ny)

      close (151)

! Save a copy that can be used as the initial condition for the
!   start of the second data assimilation loop.
      if (iter.eq.max_iter) then
        write(niter,'(i5.5)') 0
        write(iloop,'(i1.1)') 2
        fname = name1//niter//name2//iloop//name3

! Save the data.
        open(151,file=fname,&
     &    form='unformatted',access='direct',recl=4*nx*ny)

        write (151,rec=1) ((snow_d(i,j),i=1,nx),j=1,ny)
        write (151,rec=2) ((snow_depth(i,j),i=1,nx),j=1,ny)
        write (151,rec=3) ((canopy_int(i,j),i=1,nx),j=1,ny)
        write (151,rec=4) ((soft_snow_d(i,j),i=1,nx),j=1,ny)
        write (151,rec=5) ((ro_snow_grid(i,j),i=1,nx),j=1,ny)
        write (151,rec=6) ((swe_depth(i,j),i=1,nx),j=1,ny)
        write (151,rec=7) ((ro_soft_snow_old(i,j),i=1,nx),j=1,ny)
        write (151,rec=8) ((snow_d_init(i,j),i=1,nx),j=1,ny)
        write (151,rec=9) ((swe_depth_old(i,j),i=1,nx),j=1,ny)
        write (151,rec=10) ((canopy_int_old(i,j),i=1,nx),j=1,ny)
        write (151,rec=11) ((topo(i,j),i=1,nx),j=1,ny)
        write (151,rec=12) ((sum_sprec(i,j),i=1,nx),j=1,ny)

        close (151)
      endif

      return
      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine HRESTART_SAVE_DA(nx,ny,max_iter,corr_factor,&
     &  icorr_factor_index,nobs_dates)

      use snowmodel_inc
      implicit none

      integer iobs_num,nobs_dates,nx,ny,i,j,max_iter,iter
      integer icorr_factor_index(max_time_steps)
      real corr_factor(nx_max,ny_max,max_obs_dates+1)


! Save the correction factors for each observation date.
      open(152,file='hrestart/hrestart_corr_factor.gdat',&
     &  form='unformatted',access='direct',recl=4*nx*ny)

      do iobs_num=1,nobs_dates+1
        write(152,rec=iobs_num)&
     &    ((corr_factor(i,j,iobs_num),i=1,nx),j=1,ny)
      enddo

      close (152)

! Save the correction factor index.
      open(153,file='hrestart/hrestart_corr_factor_index.dat')

      do iter=1,max_iter
        write (153,*) iter,icorr_factor_index(iter)
      enddo

      close (153)

      return
      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine HRESTART_READ(nx,ny,snow_d,snow_depth,&
     &  canopy_int,soft_snow_d,ro_snow_grid,swe_depth,&
     &  ro_soft_snow_old,snow_d_init,swe_depth_old,&
     &  canopy_int_old,topo,sum_sprec,ihrestart_flag,&
     &  i_dataassim_loop)

      use snowmodel_inc
      implicit none

      integer i,j,nx,ny,ihrestart_flag,i_dataassim_loop,&
     &  i_dataassim_loop_tmp
      real snow_d(nx,ny)
      real snow_depth(nx,ny)
      real canopy_int(nx,ny)
      real soft_snow_d(nx,ny)
      real ro_snow_grid(nx,ny)
      real swe_depth(nx,ny)
      real ro_soft_snow_old(nx,ny)
      real snow_d_init(nx,ny)
      real swe_depth_old(nx,ny)
      real canopy_int_old(nx,ny)
      real topo(nx,ny)
      real sum_sprec(nx,ny)

      character*18 name1
      character*1 name2
      character*5 name3
      character*5 niter
      character*1 iloop
      character*30 fname

! Build the file name so it includes the interation number.
      name1 = 'hrestart/hrestart_'
      name2 = '_'
      name3 = '.gdat'

      if (i_dataassim_loop.lt.0.0) then
        i_dataassim_loop_tmp = 2
      else
        i_dataassim_loop_tmp = 1
      endif

      write(niter,'(i5.5)') ihrestart_flag
      write(iloop,'(i1.1)') i_dataassim_loop_tmp
      fname = name1//niter//name2//iloop//name3

! Save the data.
      open(152,file=fname,&
     &  form='unformatted',access='direct',recl=4*nx*ny)

      read (152,rec=1) ((snow_d(i,j),i=1,nx),j=1,ny)
      read (152,rec=2) ((snow_depth(i,j),i=1,nx),j=1,ny)
      read (152,rec=3) ((canopy_int(i,j),i=1,nx),j=1,ny)
      read (152,rec=4) ((soft_snow_d(i,j),i=1,nx),j=1,ny)
      read (152,rec=5) ((ro_snow_grid(i,j),i=1,nx),j=1,ny)
      read (152,rec=6) ((swe_depth(i,j),i=1,nx),j=1,ny)
      read (152,rec=7) ((ro_soft_snow_old(i,j),i=1,nx),j=1,ny)
      read (152,rec=8) ((snow_d_init(i,j),i=1,nx),j=1,ny)
      read (152,rec=9) ((swe_depth_old(i,j),i=1,nx),j=1,ny)
      read (152,rec=10) ((canopy_int_old(i,j),i=1,nx),j=1,ny)
      read (152,rec=11) ((topo(i,j),i=1,nx),j=1,ny)
      read (152,rec=12) ((sum_sprec(i,j),i=1,nx),j=1,ny)

      close (152)

      return
      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine HRESTART_READ_DA(nx,ny,max_iter,corr_factor,&
     &  icorr_factor_index,i_dataassim_loop)

      use snowmodel_inc
      implicit none

      integer iobs_num,nobs_dates,nx,ny,i,j,max_iter,iter,dummy,&
     &  i_dataassim_loop
      integer icorr_factor_index(max_time_steps)
      real corr_factor(nx_max,ny_max,max_obs_dates+1)

! Read the correction factors for each observation date.
      open(152,file='hrestart/hrestart_corr_factor.gdat',&
     &  form='unformatted',access='direct',recl=4*nx*ny)

      if (i_dataassim_loop.lt.0.0) nobs_dates = - i_dataassim_loop

      do iobs_num=1,nobs_dates+1
        read(152,rec=iobs_num)&
     &    ((corr_factor(i,j,iobs_num),i=1,nx),j=1,ny)
      enddo

      close (152)

! Read the correction factor index.
      open(153,file='hrestart/hrestart_corr_factor_index.dat')

      do iter=1,max_iter
        read (153,*) dummy,icorr_factor_index(iter)
      enddo

      close (153)

      return
      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine mk_ctl_files(nx,ny,deltax,deltay,xmn,ymn,dt,&
     &  print_inc,iyear_init,imonth_init,iday_init,xhour_init,&
     &  max_iter,undef,output_path_wo_assim,output_path_wi_assim,&
     &  print_micromet,micromet_output_fname,print_enbal,&
     &  enbal_output_fname,print_snowpack,snowpack_output_fname,&
     &  print_snowtran,snowtran_output_fname,Tabler_1_flag,&
     &  tabler_sfc_path_name,Tabler_2_flag,irun_data_assim,&
     &  print_var,print_outvars,print_multilayer,&
     &  multilayer_output_fname)

      use snowmodel_inc
      implicit none

      real deltax,deltay,dt,print_inc,xhour_init,print_micromet,&
     &  print_enbal,print_snowpack,print_snowtran,Tabler_1_flag,&
     &  Tabler_2_flag,undef,print_multilayer

      integer nx,ny,iyear_init,imonth_init,iday_init,max_iter,&
     &  irun_data_assim,k,kk

      double precision xmn,ymn

      character*1 print_var(n_print_vars)

      character*80 micromet_output_fname
      character*80 enbal_output_fname
      character*80 snowpack_output_fname
      character*80 snowtran_output_fname
      character*80 multilayer_output_fname
      character*80 tabler_sfc_path_name
      character*80 output_path
      character*80 output_path_wo_assim,output_path_wi_assim

      character*4 print_outvars(n_print_vars)
      character*80 description(n_print_vars)

! This is done below to avoid some array definition conflict (you
!   cannot use the "data" statement like this when the array is
!   being passed in and out of differet subroutines.
!     data print_outvars /'tair','relh','wspd','qsin','qlin',
!    &                    'qlem','albd','wdir','prec','rpre',
!    &                    'spre','smlt','ssub','roff','glmt',
!    &                    'snod','sden','swed','sspr','ssmt',
!    &                    'cldf','var1','var2','var3','var4',
!    &                    'var5','var6','var7','var8','var9'/

      data description / &
     &  'tair  0  0 air temperature (deg C)',&
     &  'relh  0  0 relative humidity (%)',&
     &  'wspd  0  0 wind speed (m/s)',&
     &  'qsin  0  0 incoming solar rad at the surface (W/m2)',&
     &  'qlin  0  0 incoming longwave rad at the surface (W/m2)',&

     &  'qlem  0  0 emitted longwave radiation (W/m2)',&
     &  'albd  0  0 albedo (0-1)',&
     &  'wdir  0  0 wind direction (0-360, true N)',&
     &  'prec  0  0 water-equivalent precipitation (m/time_step)',&
     &  'rpre  0  0 liquid precipitation, rain (m/time_step)',&

     &  'spre  0  0 solid precipitation, snowfall (m/time_step)',&
     &  'smlt  0  0 snow-water-equivalent melt (m)',&
     &  'ssub  0  0 static-surface sublimation (m)',&
     &  'roff  0  0 runoff from snowpack base (m/time_step)',&
     &  'glmt  0  0 snow-water-equivalent melt from glacier ice (m)',&

     &  'snod  0  0 snow depth (m)',&
     &  'sden  0  0 snow density (kg/m3)',&
     &  'swed  0  0 snow-water-equivalent depth (m)',&
     &  'sspr  0  0 summed snow precip during year (m)',&
     &  'ssmt  0  0 summed snow-water-equivalent melt (m)',&

     &  'cldf  0  0 cloud fraction (0-1)',&
     &  'var1  0  0 to be defined in future applications',&
     &  'var2  0  0 to be defined in future applications',&
     &  'var3  0  0 to be defined in future applications',&
     &  'var4  0  0 to be defined in future applications',&

     &  'var5  0  0 to be defined in future applications',&
     &  'var6  0  0 to be defined in future applications',&
     &  'var7  0  0 to be defined in future applications',&
     &  'var8  0  0 to be defined in future applications',&
     &  'var9  0  0 to be defined in future applications'/

      character*14 den_outvars(3)
      character*80 den_description(3)

      data den_outvars / &
     &  'dden_den_assim','sden_den_assim','snod_den_assim'/

      data den_description / &
     &  'dden  0  0 assimilated snow density difference (kg/m^3)',&
     &  'sden  0  0 snow density after the assimilation (kg/m^3)',&
     &  'snod  0  0 snow depth with assimilated snow density (m)'/

      print_outvars(1) = 'tair'
      print_outvars(2) = 'relh'
      print_outvars(3) = 'wspd'
      print_outvars(4) = 'qsin'
      print_outvars(5) = 'qlin'
      print_outvars(6) = 'qlem'
      print_outvars(7) = 'albd'
      print_outvars(8) = 'wdir'
      print_outvars(9) = 'prec'
      print_outvars(10) = 'rpre'
      print_outvars(11) = 'spre'
      print_outvars(12) = 'smlt'
      print_outvars(13) = 'ssub'
      print_outvars(14) = 'roff'
      print_outvars(15) = 'glmt'
      print_outvars(16) = 'snod'
      print_outvars(17) = 'sden'
      print_outvars(18) = 'swed'
      print_outvars(19) = 'sspr'
      print_outvars(20) = 'ssmt'
      print_outvars(21) = 'cldf'
      print_outvars(22) = 'var1'
      print_outvars(23) = 'var2'
      print_outvars(24) = 'var3'
      print_outvars(25) = 'var4'
      print_outvars(26) = 'var5'
      print_outvars(27) = 'var6'
      print_outvars(28) = 'var7'
      print_outvars(29) = 'var8'
      print_outvars(30) = 'var9'

! Create the GrADS .ctl (control) files to go with the GrADS
!   .gdat output files that were generated by this model run.
      if (print_micromet.eq.1.0) then 
        call mk_micromet_ctl(nx,ny,deltax,deltay,xmn,ymn,dt,&
     &    iyear_init,imonth_init,iday_init,xhour_init,max_iter,&
     &    undef,micromet_output_fname)
      endif

      if (print_enbal.eq.1.0) then
        call mk_enbal_ctl(nx,ny,deltax,deltay,xmn,ymn,dt,&
     &    iyear_init,imonth_init,iday_init,xhour_init,max_iter,&
     &    undef,enbal_output_fname)
      endif

      if (print_snowpack.eq.1.0) then
        call mk_snowpack_ctl(nx,ny,deltax,deltay,xmn,ymn,dt,&
     &    iyear_init,imonth_init,iday_init,xhour_init,max_iter,&
     &    undef,snowpack_output_fname)
      endif

      if (print_snowtran.eq.1.0) then
        call mk_snowtran_ctl(nx,ny,deltax,deltay,xmn,ymn,dt,&
     &    iyear_init,imonth_init,iday_init,xhour_init,max_iter,&
     &    undef,snowtran_output_fname)
      endif

      if (Tabler_1_flag.eq.1.0) then
        call mk_tabler_1_ctl(nx,ny,deltax,deltay,xmn,ymn,&
     &    undef,tabler_sfc_path_name)
      endif

      if (print_multilayer.eq.1.0) then
        call mk_multilayer_ctl(nx,ny,deltax,deltay,xmn,ymn,dt,&
     &    iyear_init,imonth_init,iday_init,xhour_init,max_iter,&
     &    undef,multilayer_output_fname)
      endif

! If you are not doing data assimilation, then you just need .ctl
!   files for the wo_assim directory.  If you are also doing the
!   data assimilation, then you also need .ctl files for the .gdat
!   data files that are in the wi_assim directory.
      do k=1,irun_data_assim+1
        if (k.eq.1) then
          output_path = output_path_wo_assim
        elseif (k.eq.2) then
          output_path = output_path_wi_assim
        endif

        if (Tabler_2_flag.eq.1.0) then
          call mk_tabler_2_ctl(nx,ny,deltax,deltay,xmn,ymn,dt,&
     &      iyear_init,imonth_init,iday_init,xhour_init,max_iter,&
     &      undef,output_path,k,print_inc)
        endif

        do kk=1,n_print_vars
          if (print_var(kk).eq.'y') then
            call mk_4char_vars_ctl(nx,ny,deltax,deltay,xmn,ymn,dt,&
     &        iyear_init,imonth_init,iday_init,xhour_init,max_iter,&
     &        undef,output_path,k,print_outvars(kk),description(kk),&
     &        print_inc)
          endif
        enddo

        if (print_multilayer.eq.2.0) then
          call mk_multilayer_2Dxy_ctl(nx,ny,deltax,deltay,xmn,ymn,&
     &      dt,iyear_init,imonth_init,iday_init,xhour_init,max_iter,&
     &      undef,output_path,k)

          call mk_multilayer_snod_ctl(nx,ny,deltax,deltay,xmn,ymn,&
     &      dt,iyear_init,imonth_init,iday_init,xhour_init,max_iter,&
     &      undef,output_path,k)
          call mk_multilayer_sden_ctl(nx,ny,deltax,deltay,xmn,ymn,&
     &      dt,iyear_init,imonth_init,iday_init,xhour_init,max_iter,&
     &      undef,output_path,k)
          call mk_multilayer_swed_ctl(nx,ny,deltax,deltay,xmn,ymn,&
     &      dt,iyear_init,imonth_init,iday_init,xhour_init,max_iter,&
     &      undef,output_path,k)
          call mk_multilayer_diam_ctl(nx,ny,deltax,deltay,xmn,ymn,&
     &      dt,iyear_init,imonth_init,iday_init,xhour_init,max_iter,&
     &      undef,output_path,k)
          call mk_multilayer_flux_ctl(nx,ny,deltax,deltay,xmn,ymn,&
     &      dt,iyear_init,imonth_init,iday_init,xhour_init,max_iter,&
     &      undef,output_path,k)
          call mk_multilayer_temp_ctl(nx,ny,deltax,deltay,xmn,ymn,&
     &      dt,iyear_init,imonth_init,iday_init,xhour_init,max_iter,&
     &      undef,output_path,k)
          call mk_multilayer_cond_ctl(nx,ny,deltax,deltay,xmn,ymn,&
     &      dt,iyear_init,imonth_init,iday_init,xhour_init,max_iter,&
     &      undef,output_path,k)
        endif

      enddo

!     if (irun_density_assim.eq.1) then
!       call mk_density_assim_sfc_ctl(nx,ny,deltax,deltay,xmn,ymn,
!    &    undef,den_outvars(1),den_description(1))
!       do kk=2,3
!         call mk_density_assim_ctl(nx,ny,deltax,deltay,xmn,ymn,
!    &      dt,iyear_init,imonth_init,iday_init,xhour_init,max_iter,
!    &      undef,output_path_wi_assim,print_inc,den_outvars(kk),
!    &      den_description(kk))
!       enddo
!     endif

      return
      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine mk_density_assim_sfc_ctl(nx,ny,deltax,deltay,xmn,ymn,&
     &  undef,den_outvars,den_description)

      implicit none

      integer nx,ny,igrads_dt
      double precision xmn,ymn,xmn_km,ymn_km
      real deltax,deltay,deltax_km,deltay_km,undef

      integer len_name,len_path,len_desc,trailing_blanks

      character*80 output_fname,filename,den_description
      character*80 output_path_tmp
      character*14 den_outvars
      character*3 cmo
      character*2 cdt

      data cmo /'jan'/

      deltax_km = deltax / 1000.0
      deltay_km = deltay / 1000.0
      xmn_km = xmn / 1000.0
      ymn_km = ymn / 1000.0

      igrads_dt = 1
      cdt = 'yr'

      filename = 'data/'//den_outvars//'.ctl'

      output_path_tmp = '^'

      len_path = 80 - trailing_blanks(output_path_tmp)
      output_fname = output_path_tmp(1:len_path)//den_outvars//'.gdat'
      len_name = 80 - trailing_blanks(output_fname)
      len_desc = 80 - trailing_blanks(den_description)

      open (71,file=filename)

      if (output_fname(1:1).eq.'/') then 
        write (71,50) output_fname(1:len_name)
      else 
        write (71,51) output_fname(1:len_name)
      endif

      write (71,52)
      write (71,53) undef
! (i,j) indexing.
      write (71,54) nx,1.0,1.0
      write (71,55) ny,1.0,1.0
! (meters,meters) indexing, with (zero,zero) origin.
      write (71,541) nx,0.0,deltax
      write (71,551) ny,0.0,deltay
! (km,km) indexing, with (zero,zero) origin.
      write (71,542) nx,0.0,deltax_km
      write (71,552) ny,0.0,deltay_km
! (meters,meters) indexing, with (xmn,ymn) origin.
      write (71,543) nx,xmn,deltax
      write (71,553) ny,ymn,deltay
! (km,km) indexing, with (xmn,ymn) origin.
      write (71,544) nx,xmn_km,deltax_km
      write (71,554) ny,ymn_km,deltay_km

      write (71,56)
      write (71,57) 1,12,1,cmo,9999,igrads_dt,cdt
      write (71,58)
      write (71,59) den_description(1:len_desc)
      write (71,60)

      close (71)

! This "a" by itself clips the trailing blanks in the a80 string.
   50 format ('DSET ',a)
   51 format ('DSET ^',a)
   52 format ('TITLE SnowModel Density Assimilation correction file')
   53 format ('UNDEF ',f10.1)

   54 format ('XDEF ',i8,' LINEAR ',2f20.8)
   55 format ('YDEF ',i8,' LINEAR ',2f20.8)
  541 format ('#XDEF ',i8,' LINEAR ',2f20.8)
  551 format ('#YDEF ',i8,' LINEAR ',2f20.8)
  542 format ('#XDEF ',i8,' LINEAR ',2f20.8)
  552 format ('#YDEF ',i8,' LINEAR ',2f20.8)
  543 format ('#XDEF ',i8,' LINEAR ',2f20.8)
  553 format ('#YDEF ',i8,' LINEAR ',2f20.8)
  544 format ('#XDEF ',i8,' LINEAR ',2f20.8)
  554 format ('#YDEF ',i8,' LINEAR ',2f20.8)

   56 format ('ZDEF         1 LINEAR 1 1')
! This i2.2 puts a zero in front of single digit numbers like 1.
   57 format ('TDEF ',i8,' LINEAR ',i2.2,'Z',i2.2,a3,i4,' ',i2,a2)
   58 format ('VARS     1')
   59 format (a)
   60 format ('ENDVARS')

      return
      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine mk_density_assim_ctl(nx,ny,deltax,deltay,xmn,ymn,&
     &  dt,iyear_init,imonth_init,iday_init,xhour_init,max_iter,&
     &  undef,output_path,print_inc,den_outvars,&
     &  den_description)

      implicit none

      integer nx,ny,iyear_init,imonth_init,iday_init,max_iter,&
     &  igrads_dt
      double precision xmn,ymn,xmn_km,ymn_km
      real deltax,deltay,deltax_km,deltay_km,xhour_init,dt,undef,&
     &  print_inc

      integer len_name,len_path,len_desc,trailing_blanks

      character*105 output_fname
      character*80 filename,den_description
      character*80 output_path
      character*86 output_path_tmp
      character*14 den_outvars
      character*3 cmo(12)
      character*2 cdt

      data cmo /'jan','feb','mar','apr','may','jun',&
     &          'jul','aug','sep','oct','nov','dec'/

      deltax_km = deltax / 1000.0
      deltay_km = deltay / 1000.0
      xmn_km = xmn / 1000.0
      ymn_km = ymn / 1000.0

! Convert the write interval from seconds to hours or a day.
      if (dt*print_inc.eq.86400.0) then
        igrads_dt = 1
        cdt = 'dy'
      elseif (dt*print_inc.eq.10800.0) then
        igrads_dt = 3
        cdt = 'hr'
      elseif (dt*print_inc.eq.3600.0) then
        igrads_dt = 1
        cdt = 'hr'
      else
        print *, 'This mk_ctl program has not been set up to deal'
        print *, 'with this combination of dt and print_inc values:'
        print *, 'dt =',dt,'   print_inc =',print_inc
        stop
      endif

      filename = 'ctl_files/wi_assim/'//den_outvars//'.ctl'

! Deal with the case with relative paths.
      if (output_path(1:1).ne.'/') then
        output_path_tmp = '../../'//output_path
      else
        output_path_tmp = output_path//'      '
      endif

      len_path = 86 - trailing_blanks(output_path_tmp)
      output_fname = output_path_tmp(1:len_path)//den_outvars//'.gdat'
      len_name = 105 - trailing_blanks(output_fname)
      len_desc = 80 - trailing_blanks(den_description)

      open (71,file=filename)

      if (output_fname(1:1).eq.'/') then 
        write (71,50) output_fname(1:len_name)
      else 
        write (71,51) output_fname(1:len_name)
      endif

      write (71,52)
      write (71,53) undef
! (i,j) indexing.
      write (71,54) nx,1.0,1.0
      write (71,55) ny,1.0,1.0
! (meters,meters) indexing, with (zero,zero) origin.
      write (71,541) nx,0.0,deltax
      write (71,551) ny,0.0,deltay
! (km,km) indexing, with (zero,zero) origin.
      write (71,542) nx,0.0,deltax_km
      write (71,552) ny,0.0,deltay_km
! (meters,meters) indexing, with (xmn,ymn) origin.
      write (71,543) nx,xmn,deltax
      write (71,553) ny,ymn,deltay
! (km,km) indexing, with (xmn,ymn) origin.
      write (71,544) nx,xmn_km,deltax_km
      write (71,554) ny,ymn_km,deltay_km

      write (71,56)
      write (71,57) max_iter,nint(xhour_init),iday_init,&
     &  cmo(imonth_init),iyear_init,igrads_dt,cdt
      write (71,58)
      write (71,59) den_description(1:len_desc)
      write (71,60)

      close (71)

! This "a" by itself clips the trailing blanks in the a80 string.
   50 format ('DSET ',a)
   51 format ('DSET ^',a)
   52 format ('TITLE SnowModel Density Assimilation output file')
   53 format ('UNDEF ',f10.1)

   54 format ('XDEF ',i8,' LINEAR ',2f20.8)
   55 format ('YDEF ',i8,' LINEAR ',2f20.8)
  541 format ('#XDEF ',i8,' LINEAR ',2f20.8)
  551 format ('#YDEF ',i8,' LINEAR ',2f20.8)
  542 format ('#XDEF ',i8,' LINEAR ',2f20.8)
  552 format ('#YDEF ',i8,' LINEAR ',2f20.8)
  543 format ('#XDEF ',i8,' LINEAR ',2f20.8)
  553 format ('#YDEF ',i8,' LINEAR ',2f20.8)
  544 format ('#XDEF ',i8,' LINEAR ',2f20.8)
  554 format ('#YDEF ',i8,' LINEAR ',2f20.8)

   56 format ('ZDEF         1 LINEAR 1 1')
! This i2.2 puts a zero in front of single digit numbers like 1.
   57 format ('TDEF ',i8,' LINEAR ',i2.2,'Z',i2.2,a3,i4,' ',i2,a2)
   58 format ('VARS     1')
   59 format (a)
   60 format ('ENDVARS')

      return
      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine mk_4char_vars_ctl(nx,ny,deltax,deltay,xmn,ymn,dt,&
     &  iyear_init,imonth_init,iday_init,xhour_init,max_iter,&
     &  undef,output_path,k,print_outvars,description,&
     &  print_inc)

      implicit none

      integer nx,ny,iyear_init,imonth_init,iday_init,max_iter,&
     &  igrads_dt,k
      double precision xmn,ymn,xmn_km,ymn_km
      real deltax,deltay,deltax_km,deltay_km,xhour_init,dt,undef,&
     &  print_inc

      integer len_name,len_path,len_desc,trailing_blanks

      character*95 output_fname
      character*80 filename,description
      character*80 output_path
      character*86 output_path_tmp
      character*4 print_outvars
      character*3 cmo(12)
      character*2 cdt

      data cmo /'jan','feb','mar','apr','may','jun',&
     &          'jul','aug','sep','oct','nov','dec'/

      deltax_km = deltax / 1000.0
      deltay_km = deltay / 1000.0
      xmn_km = xmn / 1000.0
      ymn_km = ymn / 1000.0

! Convert the write interval from seconds to hours or a day.
      if (dt*print_inc.eq.86400.0) then
        igrads_dt = 1
        cdt = 'dy'
      elseif (dt*print_inc.eq.10800.0) then
        igrads_dt = 3
        cdt = 'hr'
      elseif (dt*print_inc.eq.3600.0) then
        igrads_dt = 1
        cdt = 'hr'
      else
        print *, 'This mk_ctl program has not been set up to deal'
        print *, 'with this combination of dt and print_inc values:'
        print *, 'dt =',dt,'   print_inc =',print_inc
        stop
      endif

      if (k.eq.1) then
        filename = 'ctl_files/wo_assim/'//print_outvars//'.ctl'
      elseif (k.eq.2) then
        filename = 'ctl_files/wi_assim/'//print_outvars//'.ctl'
      endif

! Deal with the case with relative paths.
      if (output_path(1:1).ne.'/') then
        output_path_tmp = '../../'//output_path
      else
        output_path_tmp = output_path//'      '
      endif

      len_path = 86 - trailing_blanks(output_path_tmp)
      output_fname = output_path_tmp(1:len_path)//print_outvars//'.gdat'
      len_name = 95 - trailing_blanks(output_fname)
      len_desc = 80 - trailing_blanks(description)

      open (71,file=filename)

      if (output_fname(1:1).eq.'/') then 
        write (71,50) output_fname(1:len_name)
      else 
        write (71,51) output_fname(1:len_name)
      endif

      write (71,52)
      write (71,53) undef
! (i,j) indexing.
      write (71,54) nx,1.0,1.0
      write (71,55) ny,1.0,1.0
! (meters,meters) indexing, with (zero,zero) origin.
      write (71,541) nx,0.0,deltax
      write (71,551) ny,0.0,deltay
! (km,km) indexing, with (zero,zero) origin.
      write (71,542) nx,0.0,deltax_km
      write (71,552) ny,0.0,deltay_km
! (meters,meters) indexing, with (xmn,ymn) origin.
      write (71,543) nx,xmn,deltax
      write (71,553) ny,ymn,deltay
! (km,km) indexing, with (xmn,ymn) origin.
      write (71,544) nx,xmn_km,deltax_km
      write (71,554) ny,ymn_km,deltay_km

      write (71,56)
      write (71,57) int(real(max_iter)/print_inc),nint(xhour_init),&
     &  iday_init,cmo(imonth_init),iyear_init,igrads_dt,cdt
      write (71,58)
      write (71,59) description(1:len_desc)
      write (71,60)

      close (71)

! This "a" by itself clips the trailing blanks in the a80 string.
   50 format ('DSET ',a)
   51 format ('DSET ^',a)
   52 format ('TITLE SnowModel single-variable output file')
   53 format ('UNDEF ',f10.1)

   54 format ('XDEF ',i8,' LINEAR ',2f20.8)
   55 format ('YDEF ',i8,' LINEAR ',2f20.8)
  541 format ('#XDEF ',i8,' LINEAR ',2f20.8)
  551 format ('#YDEF ',i8,' LINEAR ',2f20.8)
  542 format ('#XDEF ',i8,' LINEAR ',2f20.8)
  552 format ('#YDEF ',i8,' LINEAR ',2f20.8)
  543 format ('#XDEF ',i8,' LINEAR ',2f20.8)
  553 format ('#YDEF ',i8,' LINEAR ',2f20.8)
  544 format ('#XDEF ',i8,' LINEAR ',2f20.8)
  554 format ('#YDEF ',i8,' LINEAR ',2f20.8)

   56 format ('ZDEF         1 LINEAR 1 1')
! This i2.2 puts a zero in front of single digit numbers like 1.
   57 format ('TDEF ',i8,' LINEAR ',i2.2,'Z',i2.2,a3,i4,' ',i2,a2)
   58 format ('VARS     1')
   59 format (a)
   60 format ('ENDVARS')

      return
      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine mk_micromet_ctl(nx,ny,deltax,deltay,xmn,ymn,dt,&
     &  iyear_init,imonth_init,iday_init,xhour_init,max_iter,&
     &  undef,output_fname)

      implicit none

      integer nx,ny,iyear_init,imonth_init,iday_init,max_iter,&
     &  igrads_dt
      double precision xmn,ymn,xmn_km,ymn_km
      real deltax,deltay,deltax_km,deltay_km,xhour_init,dt,undef

      integer len_name,trailing_blanks

      character*80 output_fname,filename
      character*83 output_fname_tmp
      character*3 cmo(12)
      character*2 cdt

      data cmo /'jan','feb','mar','apr','may','jun',&
     &          'jul','aug','sep','oct','nov','dec'/

      deltax_km = deltax / 1000.0
      deltay_km = deltay / 1000.0
      xmn_km = xmn / 1000.0
      ymn_km = ymn / 1000.0

! Convert dt from seconds to hours or a day.
      if (dt.eq.86400.0) then
        igrads_dt = 1
        cdt = 'dy'
      elseif (dt.eq.10800.0) then
        igrads_dt = 3
        cdt = 'hr'
      elseif (dt.eq.3600.0) then
        igrads_dt = 1
        cdt = 'hr'
      else
        print *, 'the mk_ctl program cannot deal with this dt value'
        stop
      endif

      filename = 'ctl_files/micromet.ctl'

! Deal with the case with relative paths.
      if (output_fname(1:1).ne.'/') then
        output_fname_tmp = '../'//output_fname
      else
        output_fname_tmp = output_fname//'   '
      endif

      len_name = 83 - trailing_blanks(output_fname_tmp)

      open (71,file=filename)

      if (output_fname(1:1).eq.'/') then 
        write (71,50) output_fname_tmp(1:len_name)
      else 
        write (71,51) output_fname_tmp(1:len_name)
      endif

      write (71,52)
      write (71,53) undef
! (i,j) indexing.
      write (71,54) nx,1.0,1.0
      write (71,55) ny,1.0,1.0
! (meters,meters) indexing, with (zero,zero) origin.
      write (71,541) nx,0.0,deltax
      write (71,551) ny,0.0,deltay
! (km,km) indexing, with (zero,zero) origin.
      write (71,542) nx,0.0,deltax_km
      write (71,552) ny,0.0,deltay_km
! (meters,meters) indexing, with (xmn,ymn) origin.
      write (71,543) nx,xmn,deltax
      write (71,553) ny,ymn,deltay
! (km,km) indexing, with (xmn,ymn) origin.
      write (71,544) nx,xmn_km,deltax_km
      write (71,554) ny,ymn_km,deltay_km

      write (71,56)
      write (71,57) max_iter,nint(xhour_init),iday_init,&
     &  cmo(imonth_init),iyear_init,igrads_dt,cdt
      write (71,58)
      write (71,59)
      write (71,60)
      write (71,61)
      write (71,62)
      write (71,63)
      write (71,64)
      write (71,65)
      write (71,66)
      write (71,67)
      write (71,68)

      close (71)

! This "a" by itself clips the trailing blanks in the a80 string.
   50 format ('DSET ',a)
   51 format ('DSET ^',a)
   52 format ('TITLE MicroMet output file')
   53 format ('UNDEF ',f10.1)

   54 format ('XDEF ',i8,' LINEAR ',2f20.8)
   55 format ('YDEF ',i8,' LINEAR ',2f20.8)
  541 format ('#XDEF ',i8,' LINEAR ',2f20.8)
  551 format ('#YDEF ',i8,' LINEAR ',2f20.8)
  542 format ('#XDEF ',i8,' LINEAR ',2f20.8)
  552 format ('#YDEF ',i8,' LINEAR ',2f20.8)
  543 format ('#XDEF ',i8,' LINEAR ',2f20.8)
  553 format ('#YDEF ',i8,' LINEAR ',2f20.8)
  544 format ('#XDEF ',i8,' LINEAR ',2f20.8)
  554 format ('#YDEF ',i8,' LINEAR ',2f20.8)

   56 format ('ZDEF         1 LINEAR 1 1')
! This i2.2 puts a zero in front of single digit numbers like 1.
   57 format ('TDEF ',i8,' LINEAR ',i2.2,'Z',i2.2,a3,i4,' ',i2,a2)
   58 format ('VARS     9')
   59 format ('tair  0  0 air temperature (deg C)')
   60 format ('relh  0  0 relative humidity (%)')
   61 format ('uwnd  0  0 meridional wind component (m/s)')
   62 format ('vwnd  0  0 zonal wind component (m/s)')
   63 format ('wspd  0  0 wind speed (m/s)')
   64 format ('wdir  0  0 wind direction (0-360, true N)')
   65 format ('qsin  0  0 incoming solar rad at the surface (W/m2)')
   66 format ('qlin  0  0 incoming longwave rad at the surface (W/m2)')
   67 format ('prec  0  0 precipitation (m/time_step)')
   68 format ('ENDVARS')

      return
      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine mk_enbal_ctl(nx,ny,deltax,deltay,xmn,ymn,dt,&
     &  iyear_init,imonth_init,iday_init,xhour_init,max_iter,&
     &  undef,output_fname)

      implicit none

      integer nx,ny,iyear_init,imonth_init,iday_init,max_iter,&
     &  igrads_dt
      double precision xmn,ymn,xmn_km,ymn_km
      real deltax,deltay,deltax_km,deltay_km,xhour_init,dt,undef

      integer len_name,trailing_blanks

      character*80 output_fname,filename
      character*83 output_fname_tmp
      character*3 cmo(12)
      character*2 cdt

      data cmo /'jan','feb','mar','apr','may','jun',&
     &          'jul','aug','sep','oct','nov','dec'/

      deltax_km = deltax / 1000.0
      deltay_km = deltay / 1000.0
      xmn_km = xmn / 1000.0
      ymn_km = ymn / 1000.0

! Convert dt from seconds to hours or a day.
      if (dt.eq.86400.0) then
        igrads_dt = 1
        cdt = 'dy'
      elseif (dt.eq.10800.0) then
        igrads_dt = 3
        cdt = 'hr'
      elseif (dt.eq.3600.0) then
        igrads_dt = 1
        cdt = 'hr'
      else
        print *, 'the mk_ctl program cannot deal with this dt value'
        stop
      endif

      filename = 'ctl_files/enbal.ctl'

! Deal with the case with relative paths.
      if (output_fname(1:1).ne.'/') then
        output_fname_tmp = '../'//output_fname
      else
        output_fname_tmp = output_fname//'   '
      endif

      len_name = 83 - trailing_blanks(output_fname_tmp)

      open (71,file=filename)

      if (output_fname(1:1).eq.'/') then 
        write (71,50) output_fname_tmp(1:len_name)
      else 
        write (71,51) output_fname_tmp(1:len_name)
      endif

      write (71,52)
      write (71,53) undef
! (i,j) indexing.
      write (71,54) nx,1.0,1.0
      write (71,55) ny,1.0,1.0
! (meters,meters) indexing, with (zero,zero) origin.
      write (71,541) nx,0.0,deltax
      write (71,551) ny,0.0,deltay
! (km,km) indexing, with (zero,zero) origin.
      write (71,542) nx,0.0,deltax_km
      write (71,552) ny,0.0,deltay_km
! (meters,meters) indexing, with (xmn,ymn) origin.
      write (71,543) nx,xmn,deltax
      write (71,553) ny,ymn,deltay
! (km,km) indexing, with (xmn,ymn) origin.
      write (71,544) nx,xmn_km,deltax_km
      write (71,554) ny,ymn_km,deltay_km

      write (71,56)
      write (71,57) max_iter,nint(xhour_init),iday_init,&
     &  cmo(imonth_init),iyear_init,igrads_dt,cdt
      write (71,58)
      write (71,59)
      write (71,60)
      write (71,61)
      write (71,62)
      write (71,63)
      write (71,64)
      write (71,65)
      write (71,66)
      write (71,67)
      write (71,68)
      write (71,69)
      write (71,70)

      close (71)

! This "a" by itself clips the trailing blanks in the a80 string.
   50 format ('DSET ',a)
   51 format ('DSET ^',a)
   52 format ('TITLE EnBal output file')
   53 format ('UNDEF ',f10.1)

   54 format ('XDEF ',i8,' LINEAR ',2f20.8)
   55 format ('YDEF ',i8,' LINEAR ',2f20.8)
  541 format ('#XDEF ',i8,' LINEAR ',2f20.8)
  551 format ('#YDEF ',i8,' LINEAR ',2f20.8)
  542 format ('#XDEF ',i8,' LINEAR ',2f20.8)
  552 format ('#YDEF ',i8,' LINEAR ',2f20.8)
  543 format ('#XDEF ',i8,' LINEAR ',2f20.8)
  553 format ('#YDEF ',i8,' LINEAR ',2f20.8)
  544 format ('#XDEF ',i8,' LINEAR ',2f20.8)
  554 format ('#YDEF ',i8,' LINEAR ',2f20.8)

   56 format ('ZDEF         1 LINEAR 1 1')
! This i2.2 puts a zero in front of single digit numbers like 1.
   57 format ('TDEF ',i8,' LINEAR ',i2.2,'Z',i2.2,a3,i4,' ',i2,a2)
   58 format ('VARS    11')
   59 format ('tair  0  0 air temperature (deg C)')
   60 format ('tsfc  0  0 surface (skin) temperature (deg C)')
   61 format ('qsin  0  0 incoming solar rad at the surface (W/m2)')
   62 format ('qlin  0  0 incoming longwave rad at the surface (W/m2)')
   63 format ('qlem  0  0 emitted longwave radiation (W/m2)')
   64 format ('qh    0  0 sensible heat flux (W/m2)')
   65 format ('qe    0  0 latent heat flux (W/m2)')
   66 format ('qc    0  0 conductive heat flux (W/m2)')
   67 format ('qm    0  0 melt energy flux (W/m2)')
   68 format ('albd  0  0 albedo (0-1)')
   69 format ('ebal  0  0 energy balance error (W/m2)')
   70 format ('ENDVARS')

      return
      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine mk_snowpack_ctl(nx,ny,deltax,deltay,xmn,ymn,dt,&
     &  iyear_init,imonth_init,iday_init,xhour_init,max_iter,&
     &  undef,output_fname)

      implicit none

      integer nx,ny,iyear_init,imonth_init,iday_init,max_iter,&
     &  igrads_dt
      double precision xmn,ymn,xmn_km,ymn_km
      real deltax,deltay,deltax_km,deltay_km,xhour_init,dt,undef

      integer len_name,trailing_blanks

      character*80 output_fname,filename
      character*83 output_fname_tmp
      character*3 cmo(12)
      character*2 cdt

      data cmo /'jan','feb','mar','apr','may','jun',&
     &          'jul','aug','sep','oct','nov','dec'/

      deltax_km = deltax / 1000.0
      deltay_km = deltay / 1000.0
      xmn_km = xmn / 1000.0
      ymn_km = ymn / 1000.0

! Convert dt from seconds to hours or a day.
      if (dt.eq.86400.0) then
        igrads_dt = 1
        cdt = 'dy'
      elseif (dt.eq.10800.0) then
        igrads_dt = 3
        cdt = 'hr'
      elseif (dt.eq.3600.0) then
        igrads_dt = 1
        cdt = 'hr'
      else
        print *, 'the mk_ctl program cannot deal with this dt value'
        stop
      endif

      filename = 'ctl_files/snowpack.ctl'

! Deal with the case with relative paths.
      if (output_fname(1:1).ne.'/') then
        output_fname_tmp = '../'//output_fname
      else
        output_fname_tmp = output_fname//'   '
      endif

      len_name = 83 - trailing_blanks(output_fname_tmp)

      open (71,file=filename)

      if (output_fname(1:1).eq.'/') then 
        write (71,50) output_fname_tmp(1:len_name)
      else 
        write (71,51) output_fname_tmp(1:len_name)
      endif

      write (71,52)
      write (71,53) undef
! (i,j) indexing.
      write (71,54) nx,1.0,1.0
      write (71,55) ny,1.0,1.0
! (meters,meters) indexing, with (zero,zero) origin.
      write (71,541) nx,0.0,deltax
      write (71,551) ny,0.0,deltay
! (km,km) indexing, with (zero,zero) origin.
      write (71,542) nx,0.0,deltax_km
      write (71,552) ny,0.0,deltay_km
! (meters,meters) indexing, with (xmn,ymn) origin.
      write (71,543) nx,xmn,deltax
      write (71,553) ny,ymn,deltay
! (km,km) indexing, with (xmn,ymn) origin.
      write (71,544) nx,xmn_km,deltax_km
      write (71,554) ny,ymn_km,deltay_km

      write (71,56)
      write (71,57) max_iter,nint(xhour_init),iday_init,&
     &  cmo(imonth_init),iyear_init,igrads_dt,cdt
      write (71,58)
      write (71,59)
      write (71,60)
      write (71,61)
      write (71,62)
      write (71,63)
      write (71,64)
      write (71,65)
      write (71,66)
      write (71,67)
      write (71,68)
      write (71,69)
      write (71,70)
      write (71,71)
      write (71,72)
      write (71,73)
      write (71,74)
      write (71,75)

      close (71)

! This "a" by itself clips the trailing blanks in the a80 string.
   50 format ('DSET ',a)
   51 format ('DSET ^',a)
   52 format ('TITLE SnowPack output file')
   53 format ('UNDEF ',f10.1)

   54 format ('XDEF ',i8,' LINEAR ',2f20.8)
   55 format ('YDEF ',i8,' LINEAR ',2f20.8)
  541 format ('#XDEF ',i8,' LINEAR ',2f20.8)
  551 format ('#YDEF ',i8,' LINEAR ',2f20.8)
  542 format ('#XDEF ',i8,' LINEAR ',2f20.8)
  552 format ('#YDEF ',i8,' LINEAR ',2f20.8)
  543 format ('#XDEF ',i8,' LINEAR ',2f20.8)
  553 format ('#YDEF ',i8,' LINEAR ',2f20.8)
  544 format ('#XDEF ',i8,' LINEAR ',2f20.8)
  554 format ('#YDEF ',i8,' LINEAR ',2f20.8)

   56 format ('ZDEF         1 LINEAR 1 1')
! This i2.2 puts a zero in front of single digit numbers like 1.
   57 format ('TDEF ',i8,' LINEAR ',i2.2,'Z',i2.2,a3,i4,' ',i2,a2)
   58 format ('VARS    16')
   59 format ('snod       0  0 snow depth (m)')
   60 format ('sden       0  0 snow density (kg/m3)')
   61 format ('swed       0  0 snow-water-equivalent depth (m)')
   62 format ('roff       0  0 runoff from snowpack base (m/time_step)')
   63 format ('rain       0  0 liquid precipitation (m/time_step)')
   64 format ('spre       0  0 solid precipitation (m/time_step)')
   65 format ('qcs        0  0 canopy sublimation (m/time_step)')
   66 format ('canopy     0  0 canopy interception store (m)')
   67 format ('sumqcs     0  0 summed canopy sublim during year (m)')
   68 format ('sumprec    0  0 summed precipitation during year (m)')
   69 format ('sumsprec   0  0 summed snow precip during year (m)')
   70 format ('sumunload  0  0 summed canopy unloading during year (m)')
   71 format ('sumroff    0  0 summed runoff during the year (m)')
   72 format ('sumswemelt 0  0 summed snow-water-equivalent melt (m)')
   73 format ('sumsublim  0  0 summed static-surface sublimation (m)')
   74 format ('wbal       0  0 summed water bal error during year (m)')
   75 format ('ENDVARS')

      return
      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine mk_snowtran_ctl(nx,ny,deltax,deltay,xmn,ymn,dt,&
     &  iyear_init,imonth_init,iday_init,xhour_init,max_iter,&
     &  undef,output_fname)

      implicit none

      integer nx,ny,iyear_init,imonth_init,iday_init,max_iter,&
     &  igrads_dt
      double precision xmn,ymn,xmn_km,ymn_km
      real deltax,deltay,deltax_km,deltay_km,xhour_init,dt,undef

      integer len_name,trailing_blanks

      character*80 output_fname,filename
      character*83 output_fname_tmp
      character*3 cmo(12)
      character*2 cdt

      data cmo /'jan','feb','mar','apr','may','jun',&
     &          'jul','aug','sep','oct','nov','dec'/

      deltax_km = deltax / 1000.0
      deltay_km = deltay / 1000.0
      xmn_km = xmn / 1000.0
      ymn_km = ymn / 1000.0

! Convert dt from seconds to hours or a day.
      if (dt.eq.86400.0) then
        igrads_dt = 1
        cdt = 'dy'
      elseif (dt.eq.10800.0) then
        igrads_dt = 3
        cdt = 'hr'
      elseif (dt.eq.3600.0) then
        igrads_dt = 1
        cdt = 'hr'
      else
        print *, 'the mk_ctl program cannot deal with this dt value'
        stop
      endif

      filename = 'ctl_files/snowtran.ctl'

! Deal with the case with relative paths.
      if (output_fname(1:1).ne.'/') then
        output_fname_tmp = '../'//output_fname
      else
        output_fname_tmp = output_fname//'   '
      endif

      len_name = 83 - trailing_blanks(output_fname_tmp)

      open (71,file=filename)

      if (output_fname(1:1).eq.'/') then 
        write (71,50) output_fname_tmp(1:len_name)
      else 
        write (71,51) output_fname_tmp(1:len_name)
      endif

      write (71,52)
      write (71,53) undef
! (i,j) indexing.
      write (71,54) nx,1.0,1.0
      write (71,55) ny,1.0,1.0
! (meters,meters) indexing, with (zero,zero) origin.
      write (71,541) nx,0.0,deltax
      write (71,551) ny,0.0,deltay
! (km,km) indexing, with (zero,zero) origin.
      write (71,542) nx,0.0,deltax_km
      write (71,552) ny,0.0,deltay_km
! (meters,meters) indexing, with (xmn,ymn) origin.
      write (71,543) nx,xmn,deltax
      write (71,553) ny,ymn,deltay
! (km,km) indexing, with (xmn,ymn) origin.
      write (71,544) nx,xmn_km,deltax_km
      write (71,554) ny,ymn_km,deltay_km

      write (71,56)
      write (71,57) max_iter,nint(xhour_init),iday_init,&
     &  cmo(imonth_init),iyear_init,igrads_dt,cdt
      write (71,58)
      write (71,59)
      write (71,60)
      write (71,61)
      write (71,62)
      write (71,63)
      write (71,64)
      write (71,65)
      write (71,66)

      close (71)

! This "a" by itself clips the trailing blanks in the a80 string.
   50 format ('DSET ',a)
   51 format ('DSET ^',a)
   52 format ('TITLE SnowTran output file')
   53 format ('UNDEF ',f10.1)

   54 format ('XDEF ',i8,' LINEAR ',2f20.8)
   55 format ('YDEF ',i8,' LINEAR ',2f20.8)
  541 format ('#XDEF ',i8,' LINEAR ',2f20.8)
  551 format ('#YDEF ',i8,' LINEAR ',2f20.8)
  542 format ('#XDEF ',i8,' LINEAR ',2f20.8)
  552 format ('#YDEF ',i8,' LINEAR ',2f20.8)
  543 format ('#XDEF ',i8,' LINEAR ',2f20.8)
  553 format ('#YDEF ',i8,' LINEAR ',2f20.8)
  544 format ('#XDEF ',i8,' LINEAR ',2f20.8)
  554 format ('#YDEF ',i8,' LINEAR ',2f20.8)

   56 format ('ZDEF         1 LINEAR 1 1')
! This i2.2 puts a zero in front of single digit numbers like 1.
   57 format ('TDEF ',i8,' LINEAR ',i2.2,'Z',i2.2,a3,i4,' ',i2,a2)
   58 format ('VARS     7')
   59 format ('snod    0  0 snow depth (m)')
   60 format ('subl    0  0 sublimation at this time step (m)')
   61 format ('salt    0  0 saltation transport at this time step (m)')
   62 format ('susp    0  0 suspended transport at this time step (m)')
   63 format ('subgrid 0  0 tabler snow redist at this time step (m)')
   64 format ('sumsubl 0  0 summed sublimation during the year (m)')
   65 format ('sumtran 0  0 summed blowing-snow transport for year (m)')
   66 format ('ENDVARS')

      return
      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine mk_tabler_1_ctl(nx,ny,deltax,deltay,xmn,ymn,&
     &  undef,output_path)

      implicit none

      integer nx,ny,igrads_dt
      double precision xmn,ymn,xmn_km,ymn_km
      real deltax,deltay,deltax_km,deltay_km,undef

      integer len_name,len_path,trailing_blanks

      character*99 output_fname
      character*80 filename
      character*80 output_path
      character*83 output_path_tmp
      character*3 cmo
      character*2 cdt

      data cmo /'jan'/

      deltax_km = deltax / 1000.0
      deltay_km = deltay / 1000.0
      xmn_km = xmn / 1000.0
      ymn_km = ymn / 1000.0

      igrads_dt = 1
      cdt = 'yr'

      filename = 'ctl_files/tabler_sfcs.ctl'

! Deal with the case with relative paths.
      if (output_path(1:1).ne.'/') then
        output_path_tmp = '../'//output_path
      else
        output_path_tmp = output_path//'   '
      endif

      len_path = 83 - trailing_blanks(output_path_tmp)
      output_fname = output_path_tmp(1:len_path)//'tabler_sfcs.gdat'
      len_name = 99 - trailing_blanks(output_fname)

      open (71,file=filename)

      if (output_fname(1:1).eq.'/') then 
        write (71,50) output_fname(1:len_name)
      else 
        write (71,51) output_fname(1:len_name)
      endif

      write (71,52)
      write (71,53) undef
! (i,j) indexing.
      write (71,54) nx,1.0,1.0
      write (71,55) ny,1.0,1.0
! (meters,meters) indexing, with (zero,zero) origin.
      write (71,541) nx,0.0,deltax
      write (71,551) ny,0.0,deltay
! (km,km) indexing, with (zero,zero) origin.
      write (71,542) nx,0.0,deltax_km
      write (71,552) ny,0.0,deltay_km
! (meters,meters) indexing, with (xmn,ymn) origin.
      write (71,543) nx,xmn,deltax
      write (71,553) ny,ymn,deltay
! (km,km) indexing, with (xmn,ymn) origin.
      write (71,544) nx,xmn_km,deltax_km
      write (71,554) ny,ymn_km,deltay_km

      write (71,56)
      write (71,57) 1,12,1,cmo,9999,igrads_dt,cdt
      write (71,58)
      write (71,59)
      write (71,60)
      write (71,61)
      write (71,62)
      write (71,63)
      write (71,64)
      write (71,65)
      write (71,66)
      write (71,67)
      write (71,68)

      close (71)

! This "a" by itself clips the trailing blanks in the a80 string.
   50 format ('DSET ',a)
   51 format ('DSET ^',a)
   52 format ('TITLE Tabler Equilibrium Surfaces for snow-free land')
   53 format ('UNDEF ',f10.1)

   54 format ('XDEF ',i8,' LINEAR ',2f20.8)
   55 format ('YDEF ',i8,' LINEAR ',2f20.8)
  541 format ('#XDEF ',i8,' LINEAR ',2f20.8)
  551 format ('#YDEF ',i8,' LINEAR ',2f20.8)
  542 format ('#XDEF ',i8,' LINEAR ',2f20.8)
  552 format ('#YDEF ',i8,' LINEAR ',2f20.8)
  543 format ('#XDEF ',i8,' LINEAR ',2f20.8)
  553 format ('#YDEF ',i8,' LINEAR ',2f20.8)
  544 format ('#XDEF ',i8,' LINEAR ',2f20.8)
  554 format ('#YDEF ',i8,' LINEAR ',2f20.8)

   56 format ('ZDEF         1 LINEAR 1 1')
! This i2.2 puts a zero in front of single digit numbers like 1.
   57 format ('TDEF ',i8,' LINEAR ',i2.2,'Z',i2.2,a3,i4,' ',i2,a2)
   58 format ('VARS     9')
   59 format ('nn    0  0  tabler surface from north winds (m)')
   60 format ('ne    0  0  tabler surface from northeast winds (m)')
   61 format ('ee    0  0  tabler surface from east winds (m)')
   62 format ('se    0  0  tabler surface from southeast winds (m)')
   63 format ('ss    0  0  tabler surface from south winds (m)')
   64 format ('sw    0  0  tabler surface from southwest winds (m)')
   65 format ('ww    0  0  tabler surface from west winds (m)')
   66 format ('nw    0  0  tabler surface from northwest winds (m)')
   67 format ('topo  0  0  topography (m)')
   68 format ('ENDVARS')

      return
      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine mk_tabler_2_ctl(nx,ny,deltax,deltay,xmn,ymn,dt,&
     &  iyear_init,imonth_init,iday_init,xhour_init,max_iter,&
     &  undef,output_path,k,print_inc)

      implicit none

      integer nx,ny,iyear_init,imonth_init,iday_init,max_iter,&
     &  igrads_dt,k
      double precision xmn,ymn,xmn_km,ymn_km
      real deltax,deltay,deltax_km,deltay_km,xhour_init,dt,undef,&
     &  print_inc

      integer len_name,len_path,trailing_blanks

      character*97 output_fname
      character*80 filename
      character*80 output_path
      character*86 output_path_tmp
      character*3 cmo(12)
      character*2 cdt

      data cmo /'jan','feb','mar','apr','may','jun',&
     &          'jul','aug','sep','oct','nov','dec'/

      deltax_km = deltax / 1000.0
      deltay_km = deltay / 1000.0
      xmn_km = xmn / 1000.0
      ymn_km = ymn / 1000.0

! Convert the write interval from seconds to hours or a day.
      if (dt*print_inc.eq.86400.0) then
        igrads_dt = 1
        cdt = 'dy'
      elseif (dt*print_inc.eq.10800.0) then
        igrads_dt = 3
        cdt = 'hr'
      elseif (dt*print_inc.eq.3600.0) then
        igrads_dt = 1
        cdt = 'hr'
      else
        print *, 'This mk_ctl program has not been set up to deal'
        print *, 'with this combination of dt and print_inc values:'
        print *, 'dt =',dt,'   print_inc =',print_inc
        stop
      endif

      if (k.eq.1) then
        filename = 'ctl_files/wo_assim/tabler_sfc_iter.ctl'
      elseif (k.eq.2) then
        filename = 'ctl_files/wi_assim/tabler_sfc_iter.ctl'
      endif

! Deal with the case with relative paths.
      if (output_path(1:1).ne.'/') then
        output_path_tmp = '../../'//output_path
      else
        output_path_tmp = output_path//'      '
      endif

      len_path = 86 - trailing_blanks(output_path_tmp)
      output_fname=output_path_tmp(1:len_path)//'tabler_sfcs_iter.gdat'
      len_name = 97 - trailing_blanks(output_fname)

      open (71,file=filename)

      if (output_fname(1:1).eq.'/') then 
        write (71,50) output_fname(1:len_name)
      else 
        write (71,51) output_fname(1:len_name)
      endif


      write (71,52)
      write (71,53) undef
! (i,j) indexing.
      write (71,54) nx,1.0,1.0
      write (71,55) ny,1.0,1.0
! (meters,meters) indexing, with (zero,zero) origin.
      write (71,541) nx,0.0,deltax
      write (71,551) ny,0.0,deltay
! (km,km) indexing, with (zero,zero) origin.
      write (71,542) nx,0.0,deltax_km
      write (71,552) ny,0.0,deltay_km
! (meters,meters) indexing, with (xmn,ymn) origin.
      write (71,543) nx,xmn,deltax
      write (71,553) ny,ymn,deltay
! (km,km) indexing, with (xmn,ymn) origin.
      write (71,544) nx,xmn_km,deltax_km
      write (71,554) ny,ymn_km,deltay_km

      write (71,56)
      write (71,57) max_iter,nint(xhour_init),iday_init,&
     &  cmo(imonth_init),iyear_init,igrads_dt,cdt
      write (71,58)
      write (71,59)
      write (71,60)

      close (71)

! This "a" by itself clips the trailing blanks in the a80 string.
   50 format ('DSET ',a)
   51 format ('DSET ^',a)
   52 format ('TITLE Tabler Equilibrium Surface output file')
   53 format ('UNDEF ',f10.1)

   54 format ('XDEF ',i8,' LINEAR ',2f20.8)
   55 format ('YDEF ',i8,' LINEAR ',2f20.8)
  541 format ('#XDEF ',i8,' LINEAR ',2f20.8)
  551 format ('#YDEF ',i8,' LINEAR ',2f20.8)
  542 format ('#XDEF ',i8,' LINEAR ',2f20.8)
  552 format ('#YDEF ',i8,' LINEAR ',2f20.8)
  543 format ('#XDEF ',i8,' LINEAR ',2f20.8)
  553 format ('#YDEF ',i8,' LINEAR ',2f20.8)
  544 format ('#XDEF ',i8,' LINEAR ',2f20.8)
  554 format ('#YDEF ',i8,' LINEAR ',2f20.8)

   56 format ('ZDEF         1 LINEAR 1 1')
! This i2.2 puts a zero in front of single digit numbers like 1.
   57 format ('TDEF ',i8,' LINEAR ',i2.2,'Z',i2.2,a3,i4,' ',i2,a2)
   58 format ('VARS     1')
   59 format ('tablersfc  0  0 Tabler equilibrium surface (m)')
   60 format ('ENDVARS')

      return
      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine mk_multilayer_ctl(nx,ny,deltax,deltay,xmn,ymn,dt,&
     &  iyear_init,imonth_init,iday_init,xhour_init,max_iter,&
     &  undef,output_fname)

      use snowmodel_inc
      implicit none

      integer nx,ny,iyear_init,imonth_init,iday_init,max_iter,&
     &  igrads_dt
      double precision xmn,ymn,xmn_km,ymn_km
      real deltax,deltay,deltax_km,deltay_km,xhour_init,dt,undef

      integer len_name,trailing_blanks

      character*80 output_fname,filename
      character*83 output_fname_tmp
      character*3 cmo(12)
      character*2 cdt

      data cmo /'jan','feb','mar','apr','may','jun',&
     &          'jul','aug','sep','oct','nov','dec'/

      deltax_km = deltax / 1000.0
      deltay_km = deltay / 1000.0
      xmn_km = xmn / 1000.0
      ymn_km = ymn / 1000.0

! Convert dt from seconds to hours or a day.
      if (dt.eq.86400.0) then
        igrads_dt = 1
        cdt = 'dy'
      elseif (dt.eq.10800.0) then
        igrads_dt = 3
        cdt = 'hr'
      elseif (dt.eq.3600.0) then
        igrads_dt = 1
        cdt = 'hr'
      else
        print *, 'the mk_ctl program cannot deal with this dt value'
        stop
      endif

      filename = 'ctl_files/multilayer.ctl'

! Deal with the case with relative paths.
      if (output_fname(1:1).ne.'/') then
        output_fname_tmp = '../'//output_fname
      else
        output_fname_tmp = output_fname//'   '
      endif

      len_name = 83 - trailing_blanks(output_fname_tmp)

      open (71,file=filename)

      if (output_fname(1:1).eq.'/') then 
        write (71,50) output_fname_tmp(1:len_name)
      else 
        write (71,51) output_fname_tmp(1:len_name)
      endif

      write (71,52)
      write (71,53) undef
! (i,j) indexing.
      write (71,54) nx,1.0,1.0
      write (71,55) ny,1.0,1.0
! (meters,meters) indexing, with (zero,zero) origin.
      write (71,541) nx,0.0,deltax
      write (71,551) ny,0.0,deltay
! (km,km) indexing, with (zero,zero) origin.
      write (71,542) nx,0.0,deltax_km
      write (71,552) ny,0.0,deltay_km
! (meters,meters) indexing, with (xmn,ymn) origin.
      write (71,543) nx,xmn,deltax
      write (71,553) ny,ymn,deltay
! (km,km) indexing, with (xmn,ymn) origin.
      write (71,544) nx,xmn_km,deltax_km
      write (71,554) ny,ymn_km,deltay_km

      write (71,56) nz_max
      write (71,57) max_iter,nint(xhour_init),iday_init,&
     &  cmo(imonth_init),iyear_init,igrads_dt,cdt
      write (71,58)
      write (71,59)
      write (71,60)
      write (71,61)
      write (71,62)
      write (71,63) nz_max
      write (71,64) nz_max
      write (71,65) nz_max
      write (71,66) nz_max
      write (71,67)

      close (71)

! This "a" by itself clips the trailing blanks in the a80 string.
   50 format ('DSET ',a)
   51 format ('DSET ^',a)
   52 format ('TITLE Multi-Layer output file')
   53 format ('UNDEF ',f10.1)

   54 format ('XDEF ',i8,' LINEAR ',2f20.8)
   55 format ('YDEF ',i8,' LINEAR ',2f20.8)
  541 format ('#XDEF ',i8,' LINEAR ',2f20.8)
  551 format ('#YDEF ',i8,' LINEAR ',2f20.8)
  542 format ('#XDEF ',i8,' LINEAR ',2f20.8)
  552 format ('#YDEF ',i8,' LINEAR ',2f20.8)
  543 format ('#XDEF ',i8,' LINEAR ',2f20.8)
  553 format ('#YDEF ',i8,' LINEAR ',2f20.8)
  544 format ('#XDEF ',i8,' LINEAR ',2f20.8)
  554 format ('#YDEF ',i8,' LINEAR ',2f20.8)

   56 format ('ZDEF    ',i6,' LINEAR 1 1')
! This i2.2 puts a zero in front of single digit numbers like 1.
   57 format ('TDEF ',i8,' LINEAR ',i2.2,'Z',i2.2,a3,i4,' ',i2,a2)
   58 format ('VARS    8')
   59 format ('KK         1  0 number of snow layers')
   60 format ('snod       1  0 snow depth (m)')
   61 format ('sden       1  0 snow density (kg/m3)')
   62 format ('swed       1  0 snow-water-equivalent depth (m)')
   63 format ('snodz ',i6,'  0 layer-specific snow depth (m)')
   64 format ('sdenz ',i6,'  0 layer-specific snow density (kg/m3)')
   65 format ('swedz ',i6,'  0 layer-specific swe depth (m)')
   66 format ('diamz ',i6,'  0 layer-specific grain diameter (m)')
   67 format ('ENDVARS')

      return
      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine mk_multilayer_2Dxy_ctl(nx,ny,deltax,deltay,xmn,ymn,&
     &  dt,iyear_init,imonth_init,iday_init,xhour_init,max_iter,&
     &  undef,output_path,k)

      use snowmodel_inc
      implicit none

      integer nx,ny,iyear_init,imonth_init,iday_init,max_iter,&
     &  igrads_dt,len_path,k
      double precision xmn,ymn,xmn_km,ymn_km
      real deltax,deltay,deltax_km,deltay_km,xhour_init,dt,undef

      integer len_name,trailing_blanks

      character*106 output_fname
      character*80 filename
      character*80 output_path
      character*86 output_path_tmp
      character*3 cmo(12)
      character*2 cdt

      data cmo /'jan','feb','mar','apr','may','jun',&
     &          'jul','aug','sep','oct','nov','dec'/

      deltax_km = deltax / 1000.0
      deltay_km = deltay / 1000.0
      xmn_km = xmn / 1000.0
      ymn_km = ymn / 1000.0

! Convert dt from seconds to hours or a day.
      if (dt.eq.86400.0) then
        igrads_dt = 1
        cdt = 'dy'
      elseif (dt.eq.10800.0) then
        igrads_dt = 3
        cdt = 'hr'
      elseif (dt.eq.3600.0) then
        igrads_dt = 1
        cdt = 'hr'
      else
        print *, 'the mk_ctl program cannot deal with this dt value'
        stop
      endif

! Define the name and location of the .ctl file.
      if (k.eq.1) then
        filename = 'ctl_files/wo_assim/multilayer_2Dxy.ctl'
      elseif (k.eq.2) then
        filename = 'ctl_files/wi_assim/multilayer_2Dxy.ctl'
      endif

! Deal with the case with relative paths.
      if (output_path(1:1).ne.'/') then
        output_path_tmp = '../../'//output_path
      else
        output_path_tmp = output_path//'      '
      endif

      len_path = 86 - trailing_blanks(output_path_tmp)
      output_fname = output_path_tmp(1:len_path)//'multilayer_2Dxy.gdat'
      len_name = 106 - trailing_blanks(output_fname)

      open (71,file=filename)

      if (output_fname(1:1).eq.'/') then 
        write (71,50) output_fname(1:len_name)
      else 
        write (71,51) output_fname(1:len_name)
      endif

      write (71,52)
      write (71,53) undef
! (i,j) indexing.
      write (71,54) nx,1.0,1.0
      write (71,55) ny,1.0,1.0
! (meters,meters) indexing, with (zero,zero) origin.
      write (71,541) nx,0.0,deltax
      write (71,551) ny,0.0,deltay
! (km,km) indexing, with (zero,zero) origin.
      write (71,542) nx,0.0,deltax_km
      write (71,552) ny,0.0,deltay_km
! (meters,meters) indexing, with (xmn,ymn) origin.
      write (71,543) nx,xmn,deltax
      write (71,553) ny,ymn,deltay
! (km,km) indexing, with (xmn,ymn) origin.
      write (71,544) nx,xmn_km,deltax_km
      write (71,554) ny,ymn_km,deltay_km

      write (71,56) nz_max
      write (71,57) max_iter,nint(xhour_init),iday_init,&
     &  cmo(imonth_init),iyear_init,igrads_dt,cdt
      write (71,58)
      write (71,59)
      write (71,60)
      write (71,61)
      write (71,62)
      write (71,67)

      close (71)

! This "a" by itself clips the trailing blanks in the a80 string.
   50 format ('DSET ',a)
   51 format ('DSET ^',a)
   52 format ('TITLE Multi-Layer output file')
   53 format ('UNDEF ',f10.1)

   54 format ('XDEF ',i8,' LINEAR ',2f20.8)
   55 format ('YDEF ',i8,' LINEAR ',2f20.8)
  541 format ('#XDEF ',i8,' LINEAR ',2f20.8)
  551 format ('#YDEF ',i8,' LINEAR ',2f20.8)
  542 format ('#XDEF ',i8,' LINEAR ',2f20.8)
  552 format ('#YDEF ',i8,' LINEAR ',2f20.8)
  543 format ('#XDEF ',i8,' LINEAR ',2f20.8)
  553 format ('#YDEF ',i8,' LINEAR ',2f20.8)
  544 format ('#XDEF ',i8,' LINEAR ',2f20.8)
  554 format ('#YDEF ',i8,' LINEAR ',2f20.8)

   56 format ('ZDEF    ',i6,' LINEAR 1 1')
! This i2.2 puts a zero in front of single digit numbers like 1.
   57 format ('TDEF ',i8,' LINEAR ',i2.2,'Z',i2.2,a3,i4,' ',i2,a2)
   58 format ('VARS    4')
   59 format ('KK         1  0 number of snow layers')
   60 format ('snod       1  0 snow depth (m)')
   61 format ('sden       1  0 snow density (kg/m3)')
   62 format ('swed       1  0 snow-water-equivalent depth (m)')
   67 format ('ENDVARS')

      return
      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine mk_multilayer_snod_ctl(nx,ny,deltax,deltay,xmn,ymn,&
     &  dt,iyear_init,imonth_init,iday_init,xhour_init,max_iter,&
     &  undef,output_path,k)

      use snowmodel_inc
      implicit none

      integer nx,ny,iyear_init,imonth_init,iday_init,max_iter,&
     &  igrads_dt,len_path,k
      double precision xmn,ymn,xmn_km,ymn_km
      real deltax,deltay,deltax_km,deltay_km,xhour_init,dt,undef

      integer len_name,trailing_blanks

      character*106 output_fname
      character*80 filename
      character*80 output_path
      character*86 output_path_tmp
      character*3 cmo(12)
      character*2 cdt

      data cmo /'jan','feb','mar','apr','may','jun',&
     &          'jul','aug','sep','oct','nov','dec'/

      deltax_km = deltax / 1000.0
      deltay_km = deltay / 1000.0
      xmn_km = xmn / 1000.0
      ymn_km = ymn / 1000.0

! Convert dt from seconds to hours or a day.
      if (dt.eq.86400.0) then
        igrads_dt = 1
        cdt = 'dy'
      elseif (dt.eq.10800.0) then
        igrads_dt = 3
        cdt = 'hr'
      elseif (dt.eq.3600.0) then
        igrads_dt = 1
        cdt = 'hr'
      else
        print *, 'the mk_ctl program cannot deal with this dt value'
        stop
      endif

! Define the name and location of the .ctl file.
      if (k.eq.1) then
        filename = 'ctl_files/wo_assim/multilayer_snod.ctl'
      elseif (k.eq.2) then
        filename = 'ctl_files/wi_assim/multilayer_snod.ctl'
      endif

! Deal with the case with relative paths.
      if (output_path(1:1).ne.'/') then
        output_path_tmp = '../../'//output_path
      else
        output_path_tmp = output_path//'      '
      endif

      len_path = 86 - trailing_blanks(output_path_tmp)
      output_fname = output_path_tmp(1:len_path)//'multilayer_snod.gdat'
      len_name = 106 - trailing_blanks(output_fname)

      open (71,file=filename)

      if (output_fname(1:1).eq.'/') then 
        write (71,50) output_fname(1:len_name)
      else 
        write (71,51) output_fname(1:len_name)
      endif

      write (71,52)
      write (71,53) undef
! (i,j) indexing.
      write (71,54) nx,1.0,1.0
      write (71,55) ny,1.0,1.0
! (meters,meters) indexing, with (zero,zero) origin.
      write (71,541) nx,0.0,deltax
      write (71,551) ny,0.0,deltay
! (km,km) indexing, with (zero,zero) origin.
      write (71,542) nx,0.0,deltax_km
      write (71,552) ny,0.0,deltay_km
! (meters,meters) indexing, with (xmn,ymn) origin.
      write (71,543) nx,xmn,deltax
      write (71,553) ny,ymn,deltay
! (km,km) indexing, with (xmn,ymn) origin.
      write (71,544) nx,xmn_km,deltax_km
      write (71,554) ny,ymn_km,deltay_km

      write (71,56) nz_max
      write (71,57) max_iter,nint(xhour_init),iday_init,&
     &  cmo(imonth_init),iyear_init,igrads_dt,cdt
      write (71,58)
      write (71,63) nz_max
      write (71,67)

      close (71)

! This "a" by itself clips the trailing blanks in the a80 string.
   50 format ('DSET ',a)
   51 format ('DSET ^',a)
   52 format ('TITLE Multi-Layer output file')
   53 format ('UNDEF ',f10.1)

   54 format ('XDEF ',i8,' LINEAR ',2f20.8)
   55 format ('YDEF ',i8,' LINEAR ',2f20.8)
  541 format ('#XDEF ',i8,' LINEAR ',2f20.8)
  551 format ('#YDEF ',i8,' LINEAR ',2f20.8)
  542 format ('#XDEF ',i8,' LINEAR ',2f20.8)
  552 format ('#YDEF ',i8,' LINEAR ',2f20.8)
  543 format ('#XDEF ',i8,' LINEAR ',2f20.8)
  553 format ('#YDEF ',i8,' LINEAR ',2f20.8)
  544 format ('#XDEF ',i8,' LINEAR ',2f20.8)
  554 format ('#YDEF ',i8,' LINEAR ',2f20.8)

   56 format ('ZDEF    ',i6,' LINEAR 1 1')
! This i2.2 puts a zero in front of single digit numbers like 1.
   57 format ('TDEF ',i8,' LINEAR ',i2.2,'Z',i2.2,a3,i4,' ',i2,a2)
   58 format ('VARS    1')
   63 format ('snodz ',i6,'  0 layer-specific snow depth (m)')
   67 format ('ENDVARS')

      return
      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine mk_multilayer_sden_ctl(nx,ny,deltax,deltay,xmn,ymn,&
     &  dt,iyear_init,imonth_init,iday_init,xhour_init,max_iter,&
     &  undef,output_path,k)

      use snowmodel_inc
      implicit none

      integer nx,ny,iyear_init,imonth_init,iday_init,max_iter,&
     &  igrads_dt,len_path,k
      double precision xmn,ymn,xmn_km,ymn_km
      real deltax,deltay,deltax_km,deltay_km,xhour_init,dt,undef

      integer len_name,trailing_blanks

      character*106 output_fname
      character*80 filename
      character*80 output_path
      character*86 output_path_tmp
      character*3 cmo(12)
      character*2 cdt

      data cmo /'jan','feb','mar','apr','may','jun',&
     &          'jul','aug','sep','oct','nov','dec'/

      deltax_km = deltax / 1000.0
      deltay_km = deltay / 1000.0
      xmn_km = xmn / 1000.0
      ymn_km = ymn / 1000.0

! Convert dt from seconds to hours or a day.
      if (dt.eq.86400.0) then
        igrads_dt = 1
        cdt = 'dy'
      elseif (dt.eq.10800.0) then
        igrads_dt = 3
        cdt = 'hr'
      elseif (dt.eq.3600.0) then
        igrads_dt = 1
        cdt = 'hr'
      else
        print *, 'the mk_ctl program cannot deal with this dt value'
        stop
      endif

! Define the name and location of the .ctl file.
      if (k.eq.1) then
        filename = 'ctl_files/wo_assim/multilayer_sden.ctl'
      elseif (k.eq.2) then
        filename = 'ctl_files/wi_assim/multilayer_sden.ctl'
      endif

! Deal with the case with relative paths.
      if (output_path(1:1).ne.'/') then
        output_path_tmp = '../../'//output_path
      else
        output_path_tmp = output_path//'      '
      endif

      len_path = 86 - trailing_blanks(output_path_tmp)
      output_fname = output_path_tmp(1:len_path)//'multilayer_sden.gdat'
      len_name = 106 - trailing_blanks(output_fname)

      open (71,file=filename)

      if (output_fname(1:1).eq.'/') then 
        write (71,50) output_fname(1:len_name)
      else 
        write (71,51) output_fname(1:len_name)
      endif

      write (71,52)
      write (71,53) undef
! (i,j) indexing.
      write (71,54) nx,1.0,1.0
      write (71,55) ny,1.0,1.0
! (meters,meters) indexing, with (zero,zero) origin.
      write (71,541) nx,0.0,deltax
      write (71,551) ny,0.0,deltay
! (km,km) indexing, with (zero,zero) origin.
      write (71,542) nx,0.0,deltax_km
      write (71,552) ny,0.0,deltay_km
! (meters,meters) indexing, with (xmn,ymn) origin.
      write (71,543) nx,xmn,deltax
      write (71,553) ny,ymn,deltay
! (km,km) indexing, with (xmn,ymn) origin.
      write (71,544) nx,xmn_km,deltax_km
      write (71,554) ny,ymn_km,deltay_km

      write (71,56) nz_max
      write (71,57) max_iter,nint(xhour_init),iday_init,&
     &  cmo(imonth_init),iyear_init,igrads_dt,cdt
      write (71,58)
      write (71,64) nz_max
      write (71,67)

      close (71)

! This "a" by itself clips the trailing blanks in the a80 string.
   50 format ('DSET ',a)
   51 format ('DSET ^',a)
   52 format ('TITLE Multi-Layer output file')
   53 format ('UNDEF ',f10.1)

   54 format ('XDEF ',i8,' LINEAR ',2f20.8)
   55 format ('YDEF ',i8,' LINEAR ',2f20.8)
  541 format ('#XDEF ',i8,' LINEAR ',2f20.8)
  551 format ('#YDEF ',i8,' LINEAR ',2f20.8)
  542 format ('#XDEF ',i8,' LINEAR ',2f20.8)
  552 format ('#YDEF ',i8,' LINEAR ',2f20.8)
  543 format ('#XDEF ',i8,' LINEAR ',2f20.8)
  553 format ('#YDEF ',i8,' LINEAR ',2f20.8)
  544 format ('#XDEF ',i8,' LINEAR ',2f20.8)
  554 format ('#YDEF ',i8,' LINEAR ',2f20.8)

   56 format ('ZDEF    ',i6,' LINEAR 1 1')
! This i2.2 puts a zero in front of single digit numbers like 1.
   57 format ('TDEF ',i8,' LINEAR ',i2.2,'Z',i2.2,a3,i4,' ',i2,a2)
   58 format ('VARS    1')
   64 format ('sdenz ',i6,'  0 layer-specific snow density (kg/m3)')
   67 format ('ENDVARS')

      return
      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine mk_multilayer_swed_ctl(nx,ny,deltax,deltay,xmn,ymn,&
     &  dt,iyear_init,imonth_init,iday_init,xhour_init,max_iter,&
     &  undef,output_path,k)

      use snowmodel_inc
      implicit none

      integer nx,ny,iyear_init,imonth_init,iday_init,max_iter,&
     &  igrads_dt,len_path,k
      double precision xmn,ymn,xmn_km,ymn_km
      real deltax,deltay,deltax_km,deltay_km,xhour_init,dt,undef

      integer len_name,trailing_blanks

      character*106 output_fname
      character*80 filename
      character*80 output_path
      character*86 output_path_tmp
      character*3 cmo(12)
      character*2 cdt

      data cmo /'jan','feb','mar','apr','may','jun',&
     &          'jul','aug','sep','oct','nov','dec'/

      deltax_km = deltax / 1000.0
      deltay_km = deltay / 1000.0
      xmn_km = xmn / 1000.0
      ymn_km = ymn / 1000.0

! Convert dt from seconds to hours or a day.
      if (dt.eq.86400.0) then
        igrads_dt = 1
        cdt = 'dy'
      elseif (dt.eq.10800.0) then
        igrads_dt = 3
        cdt = 'hr'
      elseif (dt.eq.3600.0) then
        igrads_dt = 1
        cdt = 'hr'
      else
        print *, 'the mk_ctl program cannot deal with this dt value'
        stop
      endif

! Define the name and location of the .ctl file.
      if (k.eq.1) then
        filename = 'ctl_files/wo_assim/multilayer_swed.ctl'
      elseif (k.eq.2) then
        filename = 'ctl_files/wi_assim/multilayer_swed.ctl'
      endif

! Deal with the case with relative paths.
      if (output_path(1:1).ne.'/') then
        output_path_tmp = '../../'//output_path
      else
        output_path_tmp = output_path//'      '
      endif

      len_path = 86 - trailing_blanks(output_path_tmp)
      output_fname = output_path_tmp(1:len_path)//'multilayer_swed.gdat'
      len_name = 106 - trailing_blanks(output_fname)

      open (71,file=filename)

      if (output_fname(1:1).eq.'/') then 
        write (71,50) output_fname(1:len_name)
      else 
        write (71,51) output_fname(1:len_name)
      endif

      write (71,52)
      write (71,53) undef
! (i,j) indexing.
      write (71,54) nx,1.0,1.0
      write (71,55) ny,1.0,1.0
! (meters,meters) indexing, with (zero,zero) origin.
      write (71,541) nx,0.0,deltax
      write (71,551) ny,0.0,deltay
! (km,km) indexing, with (zero,zero) origin.
      write (71,542) nx,0.0,deltax_km
      write (71,552) ny,0.0,deltay_km
! (meters,meters) indexing, with (xmn,ymn) origin.
      write (71,543) nx,xmn,deltax
      write (71,553) ny,ymn,deltay
! (km,km) indexing, with (xmn,ymn) origin.
      write (71,544) nx,xmn_km,deltax_km
      write (71,554) ny,ymn_km,deltay_km

      write (71,56) nz_max
      write (71,57) max_iter,nint(xhour_init),iday_init,&
     &  cmo(imonth_init),iyear_init,igrads_dt,cdt
      write (71,58)
      write (71,65) nz_max
      write (71,67)

      close (71)

! This "a" by itself clips the trailing blanks in the a80 string.
   50 format ('DSET ',a)
   51 format ('DSET ^',a)
   52 format ('TITLE Multi-Layer output file')
   53 format ('UNDEF ',f10.1)

   54 format ('XDEF ',i8,' LINEAR ',2f20.8)
   55 format ('YDEF ',i8,' LINEAR ',2f20.8)
  541 format ('#XDEF ',i8,' LINEAR ',2f20.8)
  551 format ('#YDEF ',i8,' LINEAR ',2f20.8)
  542 format ('#XDEF ',i8,' LINEAR ',2f20.8)
  552 format ('#YDEF ',i8,' LINEAR ',2f20.8)
  543 format ('#XDEF ',i8,' LINEAR ',2f20.8)
  553 format ('#YDEF ',i8,' LINEAR ',2f20.8)
  544 format ('#XDEF ',i8,' LINEAR ',2f20.8)
  554 format ('#YDEF ',i8,' LINEAR ',2f20.8)

   56 format ('ZDEF    ',i6,' LINEAR 1 1')
! This i2.2 puts a zero in front of single digit numbers like 1.
   57 format ('TDEF ',i8,' LINEAR ',i2.2,'Z',i2.2,a3,i4,' ',i2,a2)
   58 format ('VARS    1')
   65 format ('swedz ',i6,'  0 layer-specific swe depth (m)')
   67 format ('ENDVARS')

      return
      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine mk_multilayer_diam_ctl(nx,ny,deltax,deltay,xmn,ymn,&
     &  dt,iyear_init,imonth_init,iday_init,xhour_init,max_iter,&
     &  undef,output_path,k)

      use snowmodel_inc
      implicit none

      integer nx,ny,iyear_init,imonth_init,iday_init,max_iter,&
     &  igrads_dt,len_path,k
      double precision xmn,ymn,xmn_km,ymn_km
      real deltax,deltay,deltax_km,deltay_km,xhour_init,dt,undef

      integer len_name,trailing_blanks

      character*106 output_fname
      character*80 filename
      character*80 output_path
      character*86 output_path_tmp
      character*3 cmo(12)
      character*2 cdt

      data cmo /'jan','feb','mar','apr','may','jun',&
     &          'jul','aug','sep','oct','nov','dec'/

      deltax_km = deltax / 1000.0
      deltay_km = deltay / 1000.0
      xmn_km = xmn / 1000.0
      ymn_km = ymn / 1000.0

! Convert dt from seconds to hours or a day.
      if (dt.eq.86400.0) then
        igrads_dt = 1
        cdt = 'dy'
      elseif (dt.eq.10800.0) then
        igrads_dt = 3
        cdt = 'hr'
      elseif (dt.eq.3600.0) then
        igrads_dt = 1
        cdt = 'hr'
      else
        print *, 'the mk_ctl program cannot deal with this dt value'
        stop
      endif

! Define the name and location of the .ctl file.
      if (k.eq.1) then
        filename = 'ctl_files/wo_assim/multilayer_diam.ctl'
      elseif (k.eq.2) then
        filename = 'ctl_files/wi_assim/multilayer_diam.ctl'
      endif

! Deal with the case with relative paths.
      if (output_path(1:1).ne.'/') then
        output_path_tmp = '../../'//output_path
      else
        output_path_tmp = output_path//'      '
      endif

      len_path = 86 - trailing_blanks(output_path_tmp)
      output_fname = output_path_tmp(1:len_path)//'multilayer_diam.gdat'
      len_name = 106 - trailing_blanks(output_fname)

      open (71,file=filename)

      if (output_fname(1:1).eq.'/') then 
        write (71,50) output_fname(1:len_name)
      else 
        write (71,51) output_fname(1:len_name)
      endif

      write (71,52)
      write (71,53) undef
! (i,j) indexing.
      write (71,54) nx,1.0,1.0
      write (71,55) ny,1.0,1.0
! (meters,meters) indexing, with (zero,zero) origin.
      write (71,541) nx,0.0,deltax
      write (71,551) ny,0.0,deltay
! (km,km) indexing, with (zero,zero) origin.
      write (71,542) nx,0.0,deltax_km
      write (71,552) ny,0.0,deltay_km
! (meters,meters) indexing, with (xmn,ymn) origin.
      write (71,543) nx,xmn,deltax
      write (71,553) ny,ymn,deltay
! (km,km) indexing, with (xmn,ymn) origin.
      write (71,544) nx,xmn_km,deltax_km
      write (71,554) ny,ymn_km,deltay_km

      write (71,56) nz_max
      write (71,57) max_iter,nint(xhour_init),iday_init,&
     &  cmo(imonth_init),iyear_init,igrads_dt,cdt
      write (71,58)
      write (71,66) nz_max
      write (71,67)

      close (71)

! This "a" by itself clips the trailing blanks in the a80 string.
   50 format ('DSET ',a)
   51 format ('DSET ^',a)
   52 format ('TITLE Multi-Layer output file')
   53 format ('UNDEF ',f10.1)

   54 format ('XDEF ',i8,' LINEAR ',2f20.8)
   55 format ('YDEF ',i8,' LINEAR ',2f20.8)
  541 format ('#XDEF ',i8,' LINEAR ',2f20.8)
  551 format ('#YDEF ',i8,' LINEAR ',2f20.8)
  542 format ('#XDEF ',i8,' LINEAR ',2f20.8)
  552 format ('#YDEF ',i8,' LINEAR ',2f20.8)
  543 format ('#XDEF ',i8,' LINEAR ',2f20.8)
  553 format ('#YDEF ',i8,' LINEAR ',2f20.8)
  544 format ('#XDEF ',i8,' LINEAR ',2f20.8)
  554 format ('#YDEF ',i8,' LINEAR ',2f20.8)

   56 format ('ZDEF    ',i6,' LINEAR 1 1')
! This i2.2 puts a zero in front of single digit numbers like 1.
   57 format ('TDEF ',i8,' LINEAR ',i2.2,'Z',i2.2,a3,i4,' ',i2,a2)
   58 format ('VARS    1')
   66 format ('diamz ',i6,'  0 layer-specific grain diameter (m)')
   67 format ('ENDVARS')

      return
      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine mk_multilayer_flux_ctl(nx,ny,deltax,deltay,xmn,ymn,&
     &  dt,iyear_init,imonth_init,iday_init,xhour_init,max_iter,&
     &  undef,output_path,k)

      use snowmodel_inc
      implicit none

      integer nx,ny,iyear_init,imonth_init,iday_init,max_iter,&
     &  igrads_dt,len_path,k
      double precision xmn,ymn,xmn_km,ymn_km
      real deltax,deltay,deltax_km,deltay_km,xhour_init,dt,undef

      integer len_name,trailing_blanks

      character*106 output_fname
      character*80 filename
      character*80 output_path
      character*86 output_path_tmp
      character*3 cmo(12)
      character*2 cdt

      data cmo /'jan','feb','mar','apr','may','jun',&
     &          'jul','aug','sep','oct','nov','dec'/

      deltax_km = deltax / 1000.0
      deltay_km = deltay / 1000.0
      xmn_km = xmn / 1000.0
      ymn_km = ymn / 1000.0

! Convert dt from seconds to hours or a day.
      if (dt.eq.86400.0) then
        igrads_dt = 1
        cdt = 'dy'
      elseif (dt.eq.10800.0) then
        igrads_dt = 3
        cdt = 'hr'
      elseif (dt.eq.3600.0) then
        igrads_dt = 1
        cdt = 'hr'
      else
        print *, 'the mk_ctl program cannot deal with this dt value'
        stop
      endif

! Define the name and location of the .ctl file.
      if (k.eq.1) then
        filename = 'ctl_files/wo_assim/multilayer_flux.ctl'
      elseif (k.eq.2) then
        filename = 'ctl_files/wi_assim/multilayer_flux.ctl'
      endif

! Deal with the case with relative paths.
      if (output_path(1:1).ne.'/') then
        output_path_tmp = '../../'//output_path
      else
        output_path_tmp = output_path//'      '
      endif

      len_path = 86 - trailing_blanks(output_path_tmp)
      output_fname = output_path_tmp(1:len_path)//'multilayer_flux.gdat'
      len_name = 106 - trailing_blanks(output_fname)

      open (71,file=filename)

      if (output_fname(1:1).eq.'/') then 
        write (71,50) output_fname(1:len_name)
      else 
        write (71,51) output_fname(1:len_name)
      endif

      write (71,52)
      write (71,53) undef
! (i,j) indexing.
      write (71,54) nx,1.0,1.0
      write (71,55) ny,1.0,1.0
! (meters,meters) indexing, with (zero,zero) origin.
      write (71,541) nx,0.0,deltax
      write (71,551) ny,0.0,deltay
! (km,km) indexing, with (zero,zero) origin.
      write (71,542) nx,0.0,deltax_km
      write (71,552) ny,0.0,deltay_km
! (meters,meters) indexing, with (xmn,ymn) origin.
      write (71,543) nx,xmn,deltax
      write (71,553) ny,ymn,deltay
! (km,km) indexing, with (xmn,ymn) origin.
      write (71,544) nx,xmn_km,deltax_km
      write (71,554) ny,ymn_km,deltay_km

      write (71,56) nz_max
      write (71,57) max_iter,nint(xhour_init),iday_init,&
     &  cmo(imonth_init),iyear_init,igrads_dt,cdt
      write (71,58)
      write (71,66) nz_max
      write (71,67)

      close (71)

! This "a" by itself clips the trailing blanks in the a80 string.
   50 format ('DSET ',a)
   51 format ('DSET ^',a)
   52 format ('TITLE Multi-Layer output file')
   53 format ('UNDEF ',f10.1)

   54 format ('XDEF ',i8,' LINEAR ',2f20.8)
   55 format ('YDEF ',i8,' LINEAR ',2f20.8)
  541 format ('#XDEF ',i8,' LINEAR ',2f20.8)
  551 format ('#YDEF ',i8,' LINEAR ',2f20.8)
  542 format ('#XDEF ',i8,' LINEAR ',2f20.8)
  552 format ('#YDEF ',i8,' LINEAR ',2f20.8)
  543 format ('#XDEF ',i8,' LINEAR ',2f20.8)
  553 format ('#YDEF ',i8,' LINEAR ',2f20.8)
  544 format ('#XDEF ',i8,' LINEAR ',2f20.8)
  554 format ('#YDEF ',i8,' LINEAR ',2f20.8)

   56 format ('ZDEF    ',i6,' LINEAR 1 1')
! This i2.2 puts a zero in front of single digit numbers like 1.
   57 format ('TDEF ',i8,' LINEAR ',i2.2,'Z',i2.2,a3,i4,' ',i2,a2)
   58 format ('VARS    1')
   66 format ('fluxz ',i6,'  0 layer-specific vapor flux (kg m-2 dt-1)')
   67 format ('ENDVARS')

      return
      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine mk_multilayer_temp_ctl(nx,ny,deltax,deltay,xmn,ymn,&
     &  dt,iyear_init,imonth_init,iday_init,xhour_init,max_iter,&
     &  undef,output_path,k)

      use snowmodel_inc
      implicit none

      integer nx,ny,iyear_init,imonth_init,iday_init,max_iter,&
     &  igrads_dt,len_path,k
      double precision xmn,ymn,xmn_km,ymn_km
      real deltax,deltay,deltax_km,deltay_km,xhour_init,dt,undef

      integer len_name,trailing_blanks

      character*106 output_fname
      character*80 filename
      character*80 output_path
      character*86 output_path_tmp
      character*3 cmo(12)
      character*2 cdt

      data cmo /'jan','feb','mar','apr','may','jun',&
     &          'jul','aug','sep','oct','nov','dec'/

      deltax_km = deltax / 1000.0
      deltay_km = deltay / 1000.0
      xmn_km = xmn / 1000.0
      ymn_km = ymn / 1000.0

! Convert dt from seconds to hours or a day.
      if (dt.eq.86400.0) then
        igrads_dt = 1
        cdt = 'dy'
      elseif (dt.eq.10800.0) then
        igrads_dt = 3
        cdt = 'hr'
      elseif (dt.eq.3600.0) then
        igrads_dt = 1
        cdt = 'hr'
      else
        print *, 'the mk_ctl program cannot deal with this dt value'
        stop
      endif

! Define the name and location of the .ctl file.
      if (k.eq.1) then
        filename = 'ctl_files/wo_assim/multilayer_temp.ctl'
      elseif (k.eq.2) then
        filename = 'ctl_files/wi_assim/multilayer_temp.ctl'
      endif

! Deal with the case with relative paths.
      if (output_path(1:1).ne.'/') then
        output_path_tmp = '../../'//output_path
      else
        output_path_tmp = output_path//'      '
      endif

      len_path = 86 - trailing_blanks(output_path_tmp)
      output_fname = output_path_tmp(1:len_path)//'multilayer_temp.gdat'
      len_name = 106 - trailing_blanks(output_fname)

      open (71,file=filename)

      if (output_fname(1:1).eq.'/') then 
        write (71,50) output_fname(1:len_name)
      else 
        write (71,51) output_fname(1:len_name)
      endif

      write (71,52)
      write (71,53) undef
! (i,j) indexing.
      write (71,54) nx,1.0,1.0
      write (71,55) ny,1.0,1.0
! (meters,meters) indexing, with (zero,zero) origin.
      write (71,541) nx,0.0,deltax
      write (71,551) ny,0.0,deltay
! (km,km) indexing, with (zero,zero) origin.
      write (71,542) nx,0.0,deltax_km
      write (71,552) ny,0.0,deltay_km
! (meters,meters) indexing, with (xmn,ymn) origin.
      write (71,543) nx,xmn,deltax
      write (71,553) ny,ymn,deltay
! (km,km) indexing, with (xmn,ymn) origin.
      write (71,544) nx,xmn_km,deltax_km
      write (71,554) ny,ymn_km,deltay_km

      write (71,56) nz_max
      write (71,57) max_iter,nint(xhour_init),iday_init,&
     &  cmo(imonth_init),iyear_init,igrads_dt,cdt
      write (71,58)
      write (71,66) nz_max
      write (71,67)

      close (71)

! This "a" by itself clips the trailing blanks in the a80 string.
   50 format ('DSET ',a)
   51 format ('DSET ^',a)
   52 format ('TITLE Multi-Layer output file')
   53 format ('UNDEF ',f10.1)

   54 format ('XDEF ',i8,' LINEAR ',2f20.8)
   55 format ('YDEF ',i8,' LINEAR ',2f20.8)
  541 format ('#XDEF ',i8,' LINEAR ',2f20.8)
  551 format ('#YDEF ',i8,' LINEAR ',2f20.8)
  542 format ('#XDEF ',i8,' LINEAR ',2f20.8)
  552 format ('#YDEF ',i8,' LINEAR ',2f20.8)
  543 format ('#XDEF ',i8,' LINEAR ',2f20.8)
  553 format ('#YDEF ',i8,' LINEAR ',2f20.8)
  544 format ('#XDEF ',i8,' LINEAR ',2f20.8)
  554 format ('#YDEF ',i8,' LINEAR ',2f20.8)

   56 format ('ZDEF    ',i6,' LINEAR 1 1')
! This i2.2 puts a zero in front of single digit numbers like 1.
   57 format ('TDEF ',i8,' LINEAR ',i2.2,'Z',i2.2,a3,i4,' ',i2,a2)
   58 format ('VARS    1')
   66 format ('tempz ',i6,'  0 layer-specific snow temperature (deg K)')
   67 format ('ENDVARS')

      return
      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine mk_multilayer_cond_ctl(nx,ny,deltax,deltay,xmn,ymn,&
     &  dt,iyear_init,imonth_init,iday_init,xhour_init,max_iter,&
     &  undef,output_path,k)

      use snowmodel_inc
      implicit none

      integer nx,ny,iyear_init,imonth_init,iday_init,max_iter,&
     &  igrads_dt,len_path,k
      double precision xmn,ymn,xmn_km,ymn_km
      real deltax,deltay,deltax_km,deltay_km,xhour_init,dt,undef

      integer len_name,trailing_blanks

      character*106 output_fname
      character*80 filename
      character*80 output_path
      character*86 output_path_tmp
      character*3 cmo(12)
      character*2 cdt

      data cmo /'jan','feb','mar','apr','may','jun',&
     &          'jul','aug','sep','oct','nov','dec'/

      deltax_km = deltax / 1000.0
      deltay_km = deltay / 1000.0
      xmn_km = xmn / 1000.0
      ymn_km = ymn / 1000.0

! Convert dt from seconds to hours or a day.
      if (dt.eq.86400.0) then
        igrads_dt = 1
        cdt = 'dy'
      elseif (dt.eq.10800.0) then
        igrads_dt = 3
        cdt = 'hr'
      elseif (dt.eq.3600.0) then
        igrads_dt = 1
        cdt = 'hr'
      else
        print *, 'the mk_ctl program cannot deal with this dt value'
        stop
      endif

! Define the name and location of the .ctl file.
      if (k.eq.1) then
        filename = 'ctl_files/wo_assim/multilayer_cond.ctl'
      elseif (k.eq.2) then
        filename = 'ctl_files/wi_assim/multilayer_cond.ctl'
      endif

! Deal with the case with relative paths.
      if (output_path(1:1).ne.'/') then
        output_path_tmp = '../../'//output_path
      else
        output_path_tmp = output_path//'      '
      endif

      len_path = 86 - trailing_blanks(output_path_tmp)
      output_fname = output_path_tmp(1:len_path)//'multilayer_cond.gdat'
      len_name = 106 - trailing_blanks(output_fname)

      open (71,file=filename)

      if (output_fname(1:1).eq.'/') then 
        write (71,50) output_fname(1:len_name)
      else 
        write (71,51) output_fname(1:len_name)
      endif

      write (71,52)
      write (71,53) undef
! (i,j) indexing.
      write (71,54) nx,1.0,1.0
      write (71,55) ny,1.0,1.0
! (meters,meters) indexing, with (zero,zero) origin.
      write (71,541) nx,0.0,deltax
      write (71,551) ny,0.0,deltay
! (km,km) indexing, with (zero,zero) origin.
      write (71,542) nx,0.0,deltax_km
      write (71,552) ny,0.0,deltay_km
! (meters,meters) indexing, with (xmn,ymn) origin.
      write (71,543) nx,xmn,deltax
      write (71,553) ny,ymn,deltay
! (km,km) indexing, with (xmn,ymn) origin.
      write (71,544) nx,xmn_km,deltax_km
      write (71,554) ny,ymn_km,deltay_km

      write (71,56) nz_max
      write (71,57) max_iter,nint(xhour_init),iday_init,&
     &  cmo(imonth_init),iyear_init,igrads_dt,cdt
      write (71,58)
      write (71,66) nz_max
      write (71,67)

      close (71)

! This "a" by itself clips the trailing blanks in the a80 string.
   50 format ('DSET ',a)
   51 format ('DSET ^',a)
   52 format ('TITLE Multi-Layer output file')
   53 format ('UNDEF ',f10.1)

   54 format ('XDEF ',i8,' LINEAR ',2f20.8)
   55 format ('YDEF ',i8,' LINEAR ',2f20.8)
  541 format ('#XDEF ',i8,' LINEAR ',2f20.8)
  551 format ('#YDEF ',i8,' LINEAR ',2f20.8)
  542 format ('#XDEF ',i8,' LINEAR ',2f20.8)
  552 format ('#YDEF ',i8,' LINEAR ',2f20.8)
  543 format ('#XDEF ',i8,' LINEAR ',2f20.8)
  553 format ('#YDEF ',i8,' LINEAR ',2f20.8)
  544 format ('#XDEF ',i8,' LINEAR ',2f20.8)
  554 format ('#YDEF ',i8,' LINEAR ',2f20.8)

   56 format ('ZDEF    ',i6,' LINEAR 1 1')
! This i2.2 puts a zero in front of single digit numbers like 1.
   57 format ('TDEF ',i8,' LINEAR ',i2.2,'Z',i2.2,a3,i4,' ',i2,a2)
   58 format ('VARS    1')
   66 format ('condz ',i6,'  0 layer thermal conductivity (W m-1 K-1)')
   67 format ('ENDVARS')

      return
      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

