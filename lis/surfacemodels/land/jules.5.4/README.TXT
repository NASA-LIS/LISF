The README is a guideline for setting up LIS-JULES test cases. Please READ it
carefully before running LIS-JULES. If you want to run LIS-JULES in the same way
as US Air Force 557 WW, please obtain an rose-app-run.conf from Jerry or
Shugong. If you want to run LIS-JULES with the same settings as any UK Met
Office test case, you need a rose-app-run.conf file. Currently, the name list
files in "jules_nml" are consistent with the UK Met Office Loobos test case for
offline JULES. 

[1] rose-app-run.conf: The configuration file for the UK Met Office Rose App,
which includes FORTRAN name lists for configuring UM/JULES

  Mapping from rose-app-run.conf to JULES namelist files
  [1.1 ] [namelist:jules_elevate]       no corresponding JULES name list file 
  [1.2 ] [namelist:jules_hydrology]     ==> jules_hydrology.nml
  [1.3 ] [namelist:jules_nvegparm]      ==> nveg_params.nml
  [1.4 ] [namelist:jules_pftparm]       ==> pft_params.nml
  [1.5 ] [namelist:jules_radiation]     ==> jules_radiation.nml
  [1.6 ] [namelist:jules_sea_seaice]    no corresponding JULES name list file 
  [1.7 ] [namelist:jules_snow]          ==> jules_snow.nml
  [1.8 ] [namelist:jules_soil]          ==> jules_soil.nml
  [1.9 ] [namelist:jules_surface]       ==> jules_surface.nml
  [1.10] [namelist:jules_surface_types] ==> jules_surface_types.nml
  [1.11] [namelist:jules_triffid]       ==> triffid_params.nml
  [1.12] [namelist:jules_vegetation]    ==> jules_vegetation.nml

[2] Name list files requried by the offline JULES initialization and
configuration, which are provided in "jules_nml" directory for use. 

  [2.1] timestep.nml is for setting up JULES time step. 
    ----------------------------------------------------------------------------
    &jules_time                            | Do not change 
    l_360=.false.,                         | Do not change
    l_leap=.true.,                         | Do not change
    main_run_end='1997-12-31 23:00:00',    | Set to LIS starting date and time
    main_run_start='1996-12-31 23:00:00',  | Set to LIS ending date and time 
    timestep_len=900,                      | Set to LIS JULES timestep in seconds
    /                                      | Do not change 
    &jules_spinup                          | Do not change 
    max_spinup_cycles=0,                   | Do not change 
    /                                      | Do not change 
    ----------------------------------------------------------------------------

  [2.2] drive.nml is for forcing data configuration. The name list is created
  for feeding the "init_drive" subroutine of offline JULES.
  NOTE: Two variables in the namelist, namely z1_uv_in and z1_tq_in have to be
  double checked. 
    * z1_uv_in: Constant value for the height (m) at which the wind data are
                valid for every point. This height is relative to the
                zero-plane, not the ground. z1_uv_in should be the same as
                "JULES.5.0 reference height for forcing u and v:" in LIS
                configuration file. 
    * z1_tq_in: The height (m) at which the temperature and humidity data are
                valid) should be constant for all points. The height is relative
                to the zero-plane, not the ground. z1_tq_in should be the same
                as "JULES.5.0 reference height for forcing T and q:" in LIS
                configuration file. 
    * z1_uv_in and z1_tq_in have to be consistent with specific LIS forcing dataset.   
    * If z1_uv_in and z1_tq_in both have a default value of 10 (m) if values are
      not given in drive.nml. 

  [2.3] model_grid.nml is for setting up a dummy 1-grid domain for LIS-JULES and
  feeding the offline JULES subroutine "init_model_grid". Do not change! 

  [2.4] ancillaries.nml is for providing ancillaries for the the dummy 1-grid
  domain. Do not change!

  [2.5] initial_conditions.nml is for setting JULES initial conditions in a LIS
  ccoldstart manner. It works with file='initial_conditions.dat', please edit
  contents in 'initial_conditions.dat' accordingly. 

  [2.6] urban.nml for setting JULES urban options. Please leave it empty unless
  you are guided by JULES expert for running JULES urban model. 

  [2.7] fire.nml for setting JULES fire options. Please leave it empty unless you
  are guided by JULES expert for running JULES fire model. 

  [2.8] imogen.nml for setting JULES in imogen mode. It is not supported and
  please leave it empty. 

  [2.9] output.nml for setting output variables of offline JULES. It is not
  effective for LIS-JULES but it is needed for offline JULES configuration. Please
  leave it and don't make any changes. 

Should you have any questions or want to make any changes to these JULES name
list files, please send emails to Shugong Wang (shugong.wang@nasa.gov).

You can find detailed documentations for JULES name lists from the UK Met Office
webpage: http://jules-lsm.github.io/vn4.1/namelists/contents.html

