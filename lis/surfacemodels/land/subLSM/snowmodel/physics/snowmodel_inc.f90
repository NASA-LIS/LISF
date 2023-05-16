module snowmodel_inc
  public
  ! THESE FIRST THREE PARAMETER VALUES OFTEN NEED TO BE CHANGED
  !   FOR BIG AND LONG MODEL SIMULATIONS.

  ! nx_max and ny_max define the maximum nx and ny grid points
  !   that you can have in your model run.  If you want to run a
  !   larger domain, then you have to change these numbers and
  !   recompile the code.
  
!  integer,parameter :: nx_max=1001,ny_max=1001
!  integer,parameter :: nx_max=3131,ny_max=3603
!  integer,parameter :: nx_max=10001,ny_max=10001
!  integer,parameter :: nx_max=17636,ny_max=17784
  integer,parameter :: nx_max=22001,ny_max=22001
  
  ! max_time_steps defines the maximum number of time steps that
  !   will be used in the current compliled version of the code.
  !   If you want to run a longer time domain, then you have to
  !   change this number and recompile the code.
  
  integer, parameter :: max_time_steps=8784
  
  ! nstns_max is the maximum number of stations that can be used
  !   in the data assimilation routines.
!  integer, parameter :: nstns_max=10000
  integer, parameter :: nstns_max=22894

  ! max_obs_dates is used in the data assimilation routines.  It
  !   must be greater than or equal to the number of observation
  !   dates in the entire simulation + (plus) the number of years
  !   in the simulation.  For example, for a 6-year simulation with
  !   2 observation dates in each year, you would set max_obs_dates
  !   to be = 2obs*6yrs+6yrs = 18 or greater.  For a 6-year run with
  !   4 observation dates in 2 of the years, and 0 observation dates
  !   in the other 4 years, max_obs_dates = 8obs+6yrs = 14 or
  !   greater.
  integer, parameter :: max_obs_dates=12

  ! If you are running the multi-layer snow model (even with a single
  !   layer) nz_max must be at least one greater than max_layers in
  !   snowmodel.par.  This is because the model will build a new layer
  !   with the new snowfall and then it is merged with the layer below
  !   if you only want a single snow layer.  If you are running
  !   SnowModel's original single layer model, nz_max can be 1 (but if
  !   nz_max=2 it will avoid a warning message if you are compiling
  !   the code with gfortran).
  !     parameter (nz_max=25)
  integer, parameter :: nz_max=2

  ! This is the number of print variables that are controlled by
  !   the variable list in snowmodel.par.
  integer, parameter :: n_print_vars=30

  ! nvegtypes is the number of vegetation types used in the model
  !   run.  If you change this then you have made some big changes
  !   in the codes' vegetation description.
  integer, parameter :: nvegtypes=30
  
end module snowmodel_inc
