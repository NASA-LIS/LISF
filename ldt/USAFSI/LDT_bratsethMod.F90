!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.4
!
! Copyright (c) 2022 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!
! MODULE: LDT_bratsethMod
!
! REVISION HISTORY:
! 18 Jan 2019  Eric Kemp  Initial version forked from USAF_bratsethMod.F90
!                         in LIS.
! 05 Feb 2019  Eric Kemp  Revised to account for horizontal and vertical
!                         distance correlations following Brasnett (1999).
!                         Also, removed support for correlated observation
!                         errors.  Plus, added two new QC tests based on
!                         Brasnett.  Plus, other bug fixes.
! 02 Nov 2020  Eric Kemp  Removed blacklist at request of 557WW.
! 28 Jul 2023  Eric Kemp  Expanded station ID to 31 characters.
! 24 Aug 2023  Eric Kemp  Expanded station ID to 32 characters.
!
! DESCRIPTION:
! Source code for snow depth analysis using the Bratseth objective analysis
! scheme.
!
! USEFUL REFERENCES:
! Brasnett, B, 1999:  A global analysis of snow depth for numerical weather
!   prediction.  J Appl Meteor, 38, 726-740.
! Bratseth, A M, 1986:  Statistical interpolation by means of successive
!   corrections.  Tellus, 38A, 439-447.
! Cressie, N A C, 1993:  Statistics for Spatial Data.  Revised Edition,
!   Wiley, New York, 928 pp.
! Daley, R, 1991:  Atmospheric Data Analysis.  Cambridge University Press,
!   Cambridge, UK, 457 pp.
! Gronas, S and K H Midtbo, 1986:  Operational multivariate analyzes by
!   successive corrections.  Collection of papers presented at WMO/IUGG
!   numerical weather prediction symposium.  Tokyo, 4-8 August 1986.  J
!   Meteor Soc Japan, 64A, 61-74.
! Kalnay, E, 2003:  Atmospheric Modeling, Data Assimilation and
!   Predictability.  Cambridge University Press, Cambridge, UK, 341pp.
! Lespinas, F, V Fortin, G Roy, P Rasmussen, and T Stadnyk, 2015:  Performance
!   evaluation of the Canadian Precipitation Analysis (CaPA).  J Hydrometeor,
!   16, 2045-2064.
! Lopez, P, 2013: Experimental 4D-Var assimilation of SYNOP rain gauge data at
!   ECMWF.  Mon Wea Rev, 141, 1527-1544.
! Mahfouf, J-F, B Brasnett, and S Gagnon, 2007:  A Canadian precipitation
!   analysis (CaPA) project:  Description and preliminary results.
!   Atmos-Ocean, 45, 1-17.
! Myrick, D T, and J D Horel, 2006: Verification of surface temperature
!   forecasts from the National Digital Forecast Database over the western
!   United States.  Wea Forecasting, 21, 869-892.
! Pedder, M A, 1993:  Interpolation and filtering of spatial observations
!   using successive corrections and Gaussian filters.  Mon Wea Rev, 121,
!   2889-2902.
! Ruggiero, F H, K D Sashgeyi, A E Lipton, R V Madala, and S Raman, 1999:
!   Coupled assimilation of geostationary satellite sounder data into a
!   mesoscale model using the Bratseth approach.  Mon Wea Rev, 127, 802-821.
! Sashegyi, K D, D E Harms, R V Madala, and S Raman, 1993:  Application of
!   the Bratseth scheme for the analysis of GALE data using a mesoscale
!   model.  Mon Wea Rev, 121, 2331-2350.
!------------------------------------------------------------------------------

#include "LDT_misc.h"

module LDT_bratsethMod

   ! Defaults
   implicit none
   private

   ! Class type for storing observations, interpolated background values,
   ! observation error characteristics, and executing Bratseth analysis.
   type, public :: LDT_bratseth_t
      private
      real :: back_err_var ! Background error variance
      real :: back_err_h_corr_len ! Background error horiz correlation length
      real :: back_err_v_corr_len ! Background error vert correlation length
      integer :: nobs ! Number of observations
      character*10, allocatable :: networks(:)      ! Network name
      character*32, allocatable :: platforms(:) ! Observation station ID
      real, allocatable :: obs(:) ! Observed values
      real, allocatable :: lats(:) ! Latitude of observation (deg N)
      real, allocatable :: lons(:) ! Longitude of observation (deg E)
      real, allocatable :: elevs(:) ! Elevation of observation (m MSL)
      real, allocatable :: ob_err_vars(:) ! Error variances of observations
      real, allocatable :: backs(:) ! Background values
      integer, allocatable :: qcs(:) ! Quality control flags
      character(len=80), allocatable :: reject_reasons(:)
      ! These fields are filled during the Bratseth analysis
      real, allocatable :: inv_data_dens(:) ! Inverse data densities
      real, allocatable :: sum_ob_ests(:) ! Sum of observation estimates
      real, allocatable :: anas(:) ! Analysis at observation points
      integer :: npasses
   contains
      ! Object methods
      procedure :: new => LDT_bratseth_new
      procedure :: delete => LDT_bratseth_delete
      procedure :: count_good_obs => LDT_bratseth_count_good_obs
      procedure :: append_ob => LDT_bratseth_append_ob
      procedure :: calc_inv_data_dens => LDT_bratseth_calc_inv_data_dens
      procedure :: calc_ob_anas => LDT_bratseth_calc_ob_anas
      procedure :: calc_grid_anas => LDT_bratseth_calc_grid_anas
      procedure :: run_dup_qc => LDT_bratseth_run_dup_qc
      procedure :: run_water_qc => LDT_bratseth_run_water_qc
      procedure :: run_superstat_qc => LDT_bratseth_run_superstat_qc
      procedure :: run_back_qc => LDT_bratseth_run_back_qc
      procedure :: run_skewed_back_qc => LDT_bratseth_run_skewed_back_qc
      procedure :: run_elev_qc => LDT_bratseth_run_elev_qc
      procedure :: get_lat_lon => LDT_bratseth_get_lat_lon
      procedure :: set_back => LDT_bratseth_set_back
      procedure :: reject_ob => LDT_bratseth_reject_ob
      procedure :: resort_bad_obs => LDT_bratseth_resort_bad_obs
      procedure :: get_ob_value => LDT_bratseth_get_ob_value
      procedure :: run_nosnow_qc => LDT_bratseth_run_nosnow_qc
      procedure :: run_missing_elev_qc => LDT_bratseth_run_missing_elev_qc
      procedure :: sort_obs_by_id => LDT_bratseth_sort_obs_by_id
      procedure :: print_snowdepths => LDT_bratseth_print_snowdepths
      procedure :: count_all_obs => LDT_bratseth_count_all_obs
      procedure :: mark_good_obs => LDT_bratseth_mark_good_obs
   end type LDT_bratseth_t

   ! Private constants
   real, parameter :: MISSING = -9999

   ! Quality control flags.  Should these be public?
   real, parameter :: QC_UNKNOWN = 0
   real, parameter :: QC_GOOD = 1
   real, parameter :: QC_SUSPECT = 2
   real, parameter :: QC_REJECT = 3

contains

   ! Constructor for LDT_bratseth_t object
   subroutine LDT_bratseth_new(this,maxobs,back_err_var, &
        back_err_h_corr_len, back_err_v_corr_len)

      ! Defaults
      implicit none

      ! Arguments
      class(LDT_bratseth_t),intent(inout) :: this
      integer, intent(in) :: maxobs
      real, intent(in) :: back_err_var
      real, intent(in) :: back_err_h_corr_len
      real, intent(in) :: back_err_v_corr_len

      ! Allocate memory
      allocate(this%networks(maxobs))
      allocate(this%platforms(maxobs))
      allocate(this%obs(maxobs))
      allocate(this%lats(maxobs))
      allocate(this%lons(maxobs))
      allocate(this%elevs(maxobs))
      allocate(this%ob_err_vars(maxobs))
      allocate(this%backs(maxobs))
      allocate(this%qcs(maxobs))
      allocate(this%reject_reasons(maxobs))
      allocate(this%inv_data_dens(maxobs))
      allocate(this%sum_ob_ests(maxobs))
      allocate(this%anas(maxobs))

      ! Initialize
      this%back_err_var = back_err_var
      this%back_err_h_corr_len = back_err_h_corr_len
      this%back_err_v_corr_len = back_err_v_corr_len
      this%nobs = 0
      this%networks(:) = "NULL"
      this%platforms(:) = "NULL"
      this%obs(:) = MISSING
      this%lats(:) = MISSING
      this%lons(:) = MISSING
      this%elevs(:) = MISSING
      this%ob_err_vars(:) = MISSING
      this%backs(:) = MISSING
      this%qcs(:) = QC_REJECT ! Reset this as new obs are added
      this%reject_reasons(:) = "NONE" ! Update as QC tests are made
      this%inv_data_dens(:) = MISSING
      this%sum_ob_ests(:) = MISSING
      this%anas(:) = MISSING
      this%npasses = 0
   end subroutine LDT_bratseth_new

   ! Destructor for LDT_bratseth_t object
   subroutine LDT_bratseth_delete(this)

      ! Defaults
      implicit none

      ! Arguments
      class(LDT_bratseth_t),intent(inout) :: this

      this%back_err_var = MISSING
      this%back_err_h_corr_len = MISSING
      this%back_err_v_corr_len = MISSING
      this%nobs = 0
      if (allocated(this%networks)) deallocate(this%networks)
      if (allocated(this%platforms)) deallocate(this%platforms)
      if (allocated(this%obs)) deallocate(this%obs)
      if (allocated(this%lats)) deallocate(this%lats)
      if (allocated(this%lons)) deallocate(this%lons)
      if (allocated(this%elevs)) deallocate(this%elevs)
      if (allocated(this%ob_err_vars)) deallocate(this%ob_err_vars)
      if (allocated(this%backs)) deallocate(this%backs)
      if (allocated(this%qcs)) deallocate(this%qcs)
      if (allocated(this%reject_reasons)) deallocate(this%reject_reasons)
      if (allocated(this%inv_data_dens)) deallocate(this%inv_data_dens)
      if (allocated(this%sum_ob_ests)) deallocate(this%sum_ob_ests)
      if (allocated(this%anas)) deallocate(this%anas)
      this%npasses = 0

   end subroutine LDT_bratseth_delete

   ! Get count of all obs in LDT_bratseth_t object with "good" quality
   ! control flags
   function LDT_bratseth_count_good_obs(this) result (count)
      implicit none
      class(LDT_bratseth_t), intent(in) :: this
      integer :: count
      integer :: i
      count = 0
      do i = 1, this%nobs
         if (this%qcs(i) .ne. QC_REJECT) then
            count = count + 1
         end if
      end do ! i
   end function LDT_bratseth_count_good_obs

   ! Append a single observation into a LDT_bratseth_t object.  Value of
   ! interpolated background value is optional (useful for adding
   ! "superobservations")
   subroutine LDT_bratseth_append_ob(this, network, platform, ob, &
        lat, lon, elev, ob_err_var, back)

      ! Imports
      use LDT_logmod, only : LDT_logunit

      ! Defaults
      implicit none

      ! Arguments
      class(LDT_bratseth_t),intent(inout) :: this
      character(len=10), intent(in) :: network
      character(len=32), intent(in) :: platform
      real, intent(in) :: ob
      real, intent(in) :: lat
      real, intent(in) :: lon
      real, intent(in) :: elev
      real, intent(in) :: ob_err_var
      real, optional, intent(in) :: back

      ! Local variables
      integer :: nobs

      ! Sanity check.  Since this is intended for an operational system,
      ! just print a warning and return if we see an array bounds problem.
      nobs = this%nobs
      if (nobs .eq. size(this%obs,1)) then
         write(LDT_logunit,*) &
              '[WARN], not enough memory for assigning observation data!'
         write(LDT_logunit,*) &
              '[WARN] Can only store ',nobs,' observations'
         return
      end if

      ! Assign the value.
      nobs = nobs + 1
      this%networks(nobs) = network
      this%platforms(nobs) = platform
      this%obs(nobs) = ob
      this%lats(nobs) = lat
      ! Make sure -180 to 180 convention is used
      if (lon .gt. 180) then
         this%lons(nobs) = lon - 360.
      else
         this%lons(nobs) = lon
      end if
      this%elevs(nobs) = elev
      this%ob_err_vars(nobs) = ob_err_var
      this%qcs(nobs) = QC_UNKNOWN
      this%reject_reasons(nobs) = "NONE"
      this%nobs = nobs
      if (present(back)) then
         this%backs(nobs) = back
      end if

   end subroutine LDT_bratseth_append_ob

   ! Checks if observation network is recognized as surface snow report.
   logical function LDT_bratseth_is_snow_stn(net)
      implicit none
      character(len=10), intent(in) :: net
      logical :: answer
      answer = .false.
      if (trim(net) .eq. "AMIL") answer = .true.
      if (trim(net) .eq. "CANA") answer = .true.
      if (trim(net) .eq. "CCRHS") answer = .true.
      if (trim(net) .eq. "FAA") answer = .true.
      if (trim(net) .eq. "HADS") answer = .true.
      if (trim(net) .eq. "ICAO") answer = .true.
      if (trim(net) .eq. "MSWT") answer = .true. ! From old blacklist file
      if (trim(net) .eq. "MOBL") answer = .true.
      if (trim(net) .eq. "NWSLI") answer = .true.
      if (trim(net) .eq. "WMO") answer = .true.
      if (trim(net) .eq. "SUPEROB") answer = .true.
      LDT_bratseth_is_snow_stn = answer
   end function LDT_bratseth_is_snow_stn

   ! Calculate inverse data densities for Bratseth scheme.
   ! Since this is for LDT, no parallelization is used.
   ! NOTE:  Assumes all observations have passed QC and neighboring
   ! observations have non-correlated errors.
   subroutine LDT_bratseth_calc_inv_data_dens(this, silent)

      ! Imports
      use LDT_logMod, only: LDT_logunit

      ! Defaults
      implicit none

      ! Arguments
      class(LDT_bratseth_t), intent(inout) :: this
      logical, intent(in), optional :: silent

      ! Local variables
      integer :: num_good_obs
      real :: h_dist, v_dist, b, num, denom
      logical :: verbose
      real :: t1, t2
      integer :: i,j

      ! Control prints
      verbose = .true.
      if (present(silent)) then
         if (silent) verbose = .false.
      end if

      ! Sanity check
      num_good_obs = this%count_good_obs()
      if (num_good_obs .eq. 0) return

      if (verbose) then
         write(LDT_logunit,*) &
              '[INFO] Calculating local data densities for ', &
              num_good_obs,' obs...'
      end if

      call cpu_time(t1)

      ! Initialize inverse data densities
      this%inv_data_dens(:) = 0

      ! Loop through all unique distances between two good observations
      do j = 1, num_good_obs
         do i = 1, num_good_obs

            ! Find distance between two obs, and see if they are close enough
            ! to have correlated background errors.
            ! First, horizontal distance.
            if (i .eq. j) then
               h_dist = 0
            else
               h_dist = &
                    great_circle_distance( &
                    this%lats(i),this%lons(i), &
                    this%lats(j),this%lons(j))
            end if

            ! Next, vertical distance.
            v_dist = abs(this%elevs(i) - this%elevs(j))

            ! Calculate background error covariance
            b = back_err_cov(this%back_err_var, h_dist, v_dist, &
                 this%back_err_h_corr_len, this%back_err_v_corr_len)
            num = b ! Numerator

            ! Add observation error variance if appropriate.
            if (i .eq. j) then
               num = num + this%ob_err_vars(j)
            end if

            denom = this%back_err_var + this%ob_err_vars(i)
            denom = denom*(this%back_err_var + this%ob_err_vars(j))
            denom = sqrt(denom)

            ! NOTE:  This is data density, not inverted yet
            this%inv_data_dens(j) = this%inv_data_dens(j) + (num/denom)

         end do ! i
      end do ! j

      ! Finish data density calculations, including inversion.
      do j = 1, num_good_obs
         this%inv_data_dens(j) = &
              this%inv_data_dens(j)*(this%back_err_var + this%ob_err_vars(j))
         this%inv_data_dens(j) = &
              1. / this%inv_data_dens(j)
      end do ! j

      call cpu_time(t2)
      if (verbose) then
         write(LDT_logunit,*) &
              '[INFO] Elapsed time calculating data densities is ', &
              t2 - t1, ' seconds'
      end if

   end subroutine LDT_bratseth_calc_inv_data_dens

   ! Calculates great circle distance between two lat/lon points on globe
   ! using Vincenty formula.
   ! See https://en.wikipedia.org/wiki/Great-circle_distance
   real function great_circle_distance(lat1,lon1,lat2,lon2)

      ! Defaults
      implicit none

      ! Arguments
      real, intent(in) :: lat1, lon1, lat2, lon2

      ! Local variables
      double precision :: radlat1, radlon1, radlat2, radlon2
      double precision :: pi
      double precision :: deg2rad
      double precision :: lon_abs_diff
      double precision :: central_angle
      double precision :: term1, term2, term3

      pi = 4d0*atan(1d0)
      deg2rad = pi / 180d0
      radlat1 = dble(lat1)*deg2rad
      radlon1 = dble(lon1)*deg2rad
      radlat2 = dble(lat2)*deg2rad
      radlon2 = dble(lon2)*deg2rad
      lon_abs_diff = abs(radlon1 - radlon2)
      term1 = cos(radlat2)*sin(lon_abs_diff)
      term1 = term1*term1
      term2 = (cos(radlat1)*sin(radlat2)) - &
           (sin(radlat1)*cos(radlat2)*cos(lon_abs_diff))
      term2 = term2*term2
      term3 = (sin(radlat1)*sin(radlat2)) + &
           (cos(radlat1)*cos(radlat2)*cos(lon_abs_diff))
      central_angle = atan2( sqrt(term1 + term2) , term3 )
      great_circle_distance = real(6381000d0 * central_angle)
   end function great_circle_distance

   ! Background error covariance function
   real function back_err_cov(back_err_var, h_dist, v_dist, &
        back_err_h_corr_len, back_err_v_corr_len)
      implicit none
      real, intent(in) :: back_err_var
      real, intent(in) :: h_dist ! in meters
      real, intent(in) :: v_dist ! in meters
      real, intent(in) :: back_err_h_corr_len ! in meters
      real, intent(in) :: back_err_v_corr_len ! in meters
      real :: term1, term2, term3
      term1 = back_err_var
      term2 = soar_corr_func(h_dist, back_err_h_corr_len)
      term3 = gaussian_corr_func(v_dist, back_err_v_corr_len)
      back_err_cov = term1*term2*term3
   end function back_err_cov

   ! Second-order autoregressive correlation function
   real function soar_corr_func(dist, back_err_corr_len)
      implicit none
      real, intent(in) :: dist ! in meters
      real, intent(in) :: back_err_corr_len ! in meters
      real :: inv_corr_len, term1, term2
      inv_corr_len = 1. / back_err_corr_len
      term1 = 1 + (inv_corr_len*dist)
      term2 = exp(-1*inv_corr_len*dist)
      soar_corr_func = term1*term2
   end function soar_corr_func

   ! Gaussian correlation function
   real function gaussian_corr_func(dist, back_err_corr_len)
      implicit none
      real, intent(in) :: dist ! in meters
      real, intent(in) :: back_err_corr_len
      real :: inv_corr_len
      inv_corr_len = 1. / back_err_corr_len
      gaussian_corr_func = exp(-1*dist*dist*inv_corr_len*inv_corr_len)
   end function gaussian_corr_func

   ! Run Bratseth analysis at observation points.  Multiple interations
   ! are run until convergence is reached.  Along the way, the observation
   ! estimates from each iteration are summed for later interpolation to
   ! the grid points.  Note that the analysis is also run at the observation
   ! points because in practice the analysis converges before the iterative
   ! "observation estimates" do (see Sashegyi et al 1993).
   ! NOTE:  Assumes all observations have passed QC checks, and that
   ! observation errors are uncorrelated.
   subroutine LDT_bratseth_calc_ob_anas(this, converg_thresh, silent)

      ! Imports
      use LDT_logMod, only: LDT_logunit

      ! Defaults
      implicit none

      ! Arguments
      class(LDT_bratseth_t), intent(inout) :: this
      real, intent(in) :: converg_thresh
      logical, intent(in), optional :: silent

      ! Local variables
      integer :: num_good_obs
      real, allocatable :: prev_ests(:)
      real, allocatable :: prev_anas(:)
      real, allocatable :: new_anas(:)
      real, allocatable :: new_ests(:)
      real :: h_dist, v_dist, b, weight
      real :: y_new, y_prev
      logical :: done
      integer :: icount
      integer :: i,j
      logical :: verbose
      real :: t0, t1, t2
      real :: rmsd, maxrmsd
      integer :: imaxrmsd

      verbose = .true.
      if (present(silent)) then
         if (silent) verbose = .false.
      end if

      if (verbose) then
         write(LDT_logunit,*) &
              '[INFO] Running analysis at observation points...'
      end if

      ! Sanity checks
      num_good_obs = this%count_good_obs()
      if (num_good_obs .eq. 0) return

      ! Perform analysi at observation points.  See Bratseth (1986) or
      ! Sashegyi et al (1993).
      !
      ! In each iteration, we calculate both an updated observation estimate
      ! and an updated analysis.  The previous observation estimate vector
      ! is used in both calculations, but the interpolation weights differ
      ! which cause the observation estimates and analysis values to drift
      ! apart with each iteration.  We use the change in analysis values to
      ! see if the analysis has converged; if so, we terminate the iterations.
      !
      ! The output we need are the summed observation estimates and the number
      ! of iterations (passes) required for convergence; both are used later
      ! to interpolate the analysis to the LIS grid points.
      call cpu_time(t0)
      allocate(new_anas(num_good_obs))
      allocate(new_ests(num_good_obs))
      allocate(prev_anas(num_good_obs))
      allocate(prev_ests(num_good_obs))

      prev_ests(:) = this%backs(1:num_good_obs) ! First guess
      prev_anas(:) = this%backs(1:num_good_obs) ! First guess
      this%sum_ob_ests(:) = 0
      this%npasses = 0

      do ! Iterate until convergence
         call cpu_time(t1)
         new_anas(:) = 0
         new_ests(:) = 0
         do j = 1, num_good_obs
            do i = 1, num_good_obs

               ! Calculate the distance between the obs, and see if they
               ! are close enough to have correlated background errors.
               ! First, horizontal distance
               if (i .eq. j) then
                  h_dist = 0
               else
                  h_dist = &
                       great_circle_distance( &
                       this%lats(i),this%lons(i), &
                       this%lats(j),this%lons(j))
               end if

               ! Next, vertical distance.
               v_dist = abs(this%elevs(i) - this%elevs(j))

               b = back_err_cov(this%back_err_var, h_dist, v_dist, &
                    this%back_err_h_corr_len, this%back_err_v_corr_len)

               ! First, update the observation estimate at point j
               weight = b
               if (i .eq. j) then
                  weight = weight + this%ob_err_vars(i)
               end if
               weight = weight * this%inv_data_dens(i)
               new_ests(j) = new_ests(j) + &
                    (weight*(this%obs(i) - prev_ests(i)))

               ! Second, update the analysis at the observation point j
               weight = b * this%inv_data_dens(i)
               new_anas(j) = new_anas(j) + &
                    (weight*(this%obs(i) - prev_ests(i)))

            end do ! i
         end do ! j

         ! Finish analysis and observation estimates for this iteration.
         ! Update sum of observation estimates
         do j = 1, num_good_obs
            new_ests(j) = prev_ests(j) + new_ests(j)
            new_anas(j) = prev_anas(j) + new_anas(j)
            this%sum_ob_ests(j) = this%sum_ob_ests(j) + prev_ests(j)
         end do ! j
         this%npasses = this%npasses + 1

         ! See if analysis has converged by comparing current and previous
         ! analysis values.  In practice, the analysis will converge before
         ! the observation estimates do (unless observations are assumed
         ! perfect).  Then, overwrite previous values.
         done = .true. ! We will change this below if necessary.
         maxrmsd = 0
         imaxrmsd = 0
         do j = 1, num_good_obs

            ! Check for convergence
            y_prev = prev_anas(j)
            y_new = new_anas(j)
            rmsd = sqrt((y_prev - y_new)*(y_prev - y_new))
            if (rmsd .gt. converg_thresh) then
               if (rmsd .gt. maxrmsd) then
                  maxrmsd = rmsd
                  imaxrmsd = j
               end if
               done = .false.
            endif

            ! Overwrite previous values
            prev_ests(j) = new_ests(j)
            prev_anas(j) = new_anas(j)

         end do ! j

         if (done) exit ! Analysis has converged

         if (verbose) then
            write(LDT_logunit,*) &
                 '[INFO] Bratseth scheme not converged yet after ',&
                 this%npasses, 'iterations'
            write(LDT_logunit,*) &
                 '[INFO] Max RMS change ',maxrmsd,' at i = ', imaxrmsd
            write(LDT_logunit,*) &
                 '[INFO] network, platform, obs, back, ana, est, ', &
                 'dataDensity: ', &
                 trim(this%networks(imaxrmsd)), ' ',&
                 trim(this%platforms(imaxrmsd)), &
                 ' ',this%obs(imaxrmsd),&
                 ' ',this%backs(imaxrmsd),&
                 ' ',new_anas(imaxrmsd),&
                 ' ',new_ests(imaxrmsd),&
                 ' ',1./this%inv_data_dens(imaxrmsd)
         end if ! verbose

         call cpu_time(t2)
         if (verbose) then
            write(LDT_logunit,*) &
                 '[INFO] Elapsed time for this iteration is ',t2 - t1, &
                 ' seconds'
         end if

         ! Escape if too many iterations.  Since this is intended for
         ! production, we will allow the program to continue instead of
         ! aborting.
         if (this%npasses .eq. 100) then
            write(LDT_logunit,*) &
                 '[WARN] Bratseth failed to converge after ', &
                 this%npasses, ' iterations!'
            write(LDT_logunit,*) &
                 '[WARN] Will stop iterating'
            flush(LDT_logunit)
            exit
         end if

      end do ! Iterate until convergence or too many passes

      ! Save the local analysis values
      do j = 1, num_good_obs
         this%anas(j) = new_anas(j)
      end do ! j

      if (verbose) then
         if (done) then
            write(LDT_logunit,*) &
                 '[INFO] Bratseth analysis converged after ',&
                 this%npasses, ' iterations'
         end if
      end if

      call cpu_time(t2)
      if (verbose) then
         write(LDT_logunit,*) &
              '[INFO] Elapsed time for final iteration is ',t2 - t1, &
              ' seconds'
         write(LDT_logunit,*) &
              '[INFO] Total elapsed time for all iterations is ',t2 - t0, &
              ' seconds'
         write(LDT_logunit,*) &
              '---------------------------------------------------------------'
      end if

      if (verbose) then
         icount = num_good_obs
         write(LDT_logunit,*) '[INFO] ',icount,' obs used in this analysis'
      end if

      ! Clean up
      deallocate(new_anas)
      deallocate(new_ests)
      deallocate(prev_anas)
      deallocate(prev_ests)

   end subroutine LDT_bratseth_calc_ob_anas

   ! Perform Bratseth analysis at grid points.  Assumes Bratseth scheme
   ! was already run at the observation points.
   !
   ! The interpolation from observation points to grid points is done in
   ! a single pass, similar to Daley (1991) or Kalnay (2003).  This greatly
   ! saves time compared to the original Bratseth (1986) or Sashegyi et al
   ! (1993) approaches, where ob-to-grid interpolation was done as each
   ! analysis pass was performed at the observation points.
   !
   ! NOTE:  Bratseth values are not interpolated to water points.
   ! NOTE:  We assume all observations passed QC tests.
   subroutine LDT_bratseth_calc_grid_anas(this,n,ncols,nrows,gbacks,gelevs, &
        ganas,skip_grid_points)

      ! Imports
      use LDT_coreMod, only: LDT_domain
      use LDT_logMod, only: LDT_logunit

      ! Defaults
      implicit none

      ! Arguments
      class(LDT_bratseth_t), intent(in) :: this
      integer, intent(in) :: n
      integer, intent(in) :: ncols
      integer, intent(in) :: nrows
      real, intent(in) :: gbacks(ncols,nrows)
      real, intent(in) :: gelevs(ncols,nrows)
      real, intent(inout) :: ganas(ncols,nrows)
      logical,intent(in) :: skip_grid_points(ncols,nrows)

      ! Local variables
      real :: t1, t2
      integer :: num_good_obs
      real :: h_dist, v_dist, weight
      integer :: c,r,j
      real :: tmp_ana
      real :: lat, lon

      ! Sanity check
      num_good_obs = this%count_good_obs()
      if (num_good_obs .eq. 0) return

      write(LDT_logunit,*) &
           '[INFO] Calculating analysis at grid points...'
      call cpu_time(t1)

      ! Now calculate the analysis at each grid point.
      do r = 1, nrows
         do c = 1, ncols

            if (skip_grid_points(c,r)) cycle

            tmp_ana = gbacks(c,r)

            ! Find lat/lon of tile
            lat = LDT_domain(n)%lat(c+(r-1)*ncols)
            lon = LDT_domain(n)%lon(c+(r-1)*ncols)

            do j = 1, num_good_obs

               ! See if observation is close enough for background error
               ! to be correlated.
               ! First, horizontal distance.
               h_dist = great_circle_distance(lat,lon, &
                     this%lats(j), this%lons(j))

               ! Next, vertical distance
               v_dist = abs(gelevs(c,r) - this%elevs(j))

               weight = &
                    back_err_cov(this%back_err_var, h_dist, v_dist, &
                    this%back_err_h_corr_len, this%back_err_v_corr_len)
               weight = weight * this%inv_data_dens(j)

               tmp_ana = tmp_ana + &
                    (weight * ( (this%npasses*this%obs(j)) - &
                                (this%sum_ob_ests(j)     ) ))

            end do ! j

            ganas(c,r) = tmp_ana

         end do ! c
      end do ! r

      call cpu_time(t2)
      write(LDT_logunit,*)'[INFO] Elapsed time for grid analysis is ', &
           t2 - t1,' seconds'

   end subroutine LDT_bratseth_calc_grid_anas

   ! QC checks for duplicate gage reports.  Based on Mahfouf et al (2007).
   !
   ! If duplicates are found for a particular station but all are identical,
   ! only one report is preserved and the rest are rejected.  Otherwise,
   ! if two different reports from the same station are found, a superob will
   ! be created if the spread is smaller than the observation error variance
   ! and both original reports will be rejected.  If more than two
   ! unique reports exist for the same station, all will be rejected.
   subroutine LDT_bratseth_run_dup_qc(this)

      ! Imports
      use LDT_logMod, only: LDT_logunit

      ! Defaults
      implicit none

      ! Arguments
      class(LDT_bratseth_t), intent(inout) :: this

      ! Local variables
      type close_obs ! Linked list type
         integer :: ob_index
         type(close_obs), pointer :: next
      end type close_obs
      type(close_obs), pointer :: head, tail, new, ptr
      integer :: num_good_obs
      integer :: total_reject_count
      integer :: total_merge_count
      integer :: total_create_count
      integer :: count_dups
      real :: mean, back, newlat, newlon, newelev
      real :: ob_err_var, ob_err_corr_len
      character(len=10) :: network
      character(len=32) :: platform
      real :: diff
      integer :: i,j
      logical :: reject_all
      real :: t1, t2
      logical :: location_issue

      ! Sanity check
      num_good_obs = this%count_good_obs()
      if (num_good_obs .eq. 0) then
         write(LDT_logunit,*) &
              '[INFO] dupQC found no good observations to test'
         return
      end if

      call cpu_time(t1)

      total_reject_count = 0
      total_merge_count = 0
      total_create_count = 0

      do j = 1, num_good_obs

         count_dups = 0
         mean = 0
         nullify(head,tail) ! Initialize linked list

         do i = j+1, num_good_obs

            if (this%qcs(i) .eq. QC_REJECT) cycle
            if (trim(this%networks(i)) .ne. trim(this%networks(j))) cycle
            if (trim(this%platforms(i)) .ne. trim(this%platforms(j))) cycle

            ! Duplicate found.  Store in linked list.
            count_dups = count_dups + 1
            allocate(new)
            new%ob_index = i
            nullify(new%next)
            if (associated(head)) then
               tail%next => new
               tail => new
            else ! First entry
               head => new
               tail => new
            end if

         end do ! i

         ! If we have no duplicates, just move on.
         if (count_dups .eq. 0) cycle

         ! Handle the special case of inconsistent locations.  In this event,
         ! just save the first ob.
         location_issue = .false.
         if (count_dups .gt. 0) then

            ! First pass
            ptr => head
            do i = 1, count_dups
               diff = this%lats(ptr%ob_index) - this%lats(i)
               if (diff .ne. 0) location_issue = .true.
               diff = this%lons(ptr%ob_index) - this%lons(i)
               if (diff .ne. 0) location_issue = .true.
               if (location_issue) exit ! Get out of loop
               ptr => ptr%next
            end do ! i

            if (location_issue) then
               ! Toss all duplicate reports, just save the first one.
               ptr => head
               do i = 1, count_dups
                  this%qcs(ptr%ob_index) = QC_REJECT
                  this%reject_reasons(ptr%ob_index) = "DUP_QC REJECT"
                  total_reject_count = total_reject_count + 1
                  ! write(LDT_logunit,*) &
!                        '[INFO] dupQC rejecting ob w/ changing lat/lon  r: ', &
!                        ptr%ob_index, &
!                        ' net: ',trim(this%networks(ptr%ob_index)), &
!                        ' platform: ',trim(this%platforms(ptr%ob_index)), &
!                        ' lat: ',this%lats(ptr%ob_index), &
!                        ' lon: ',this%lons(ptr%ob_index), &
!                        ' elev: ',this%elevs(ptr%ob_index), &
!                        ' obs: ',this%obs(ptr%ob_index), &
!                        ' back: ',this%backs(ptr%ob_index)
!                   write(LDT_logunit,*) &
!                        '------------------------------------------------------'
                  ptr => ptr%next
               end do ! i
            end if
         end if

         ! If we have more than one duplicate, determining a sane procedure
         ! for handling all the reports becomes difficult.  We'll keep it
         ! simple and reject all reports if more than one unique measurement
         ! amount is found; otherwise, we'll just reject the exact duplicates.
         if (count_dups .gt. 1 .and. .not. location_issue) then
            ptr => head
            reject_all = .false.
            do i = 1, count_dups
               diff = this%obs(ptr%ob_index) - this%obs(j)
               if (diff .eq. 0) then
                  this%qcs(ptr%ob_index) = QC_REJECT
                  this%reject_reasons(ptr%ob_index) = "DUP_QC REJECT"
                  ! write(LDT_logunit,*) &
!                        '[INFO] dupQC rejecting exact duplicate1 ob i: ', &
!                        ptr%ob_index, &
!                        ' network: ',trim(this%networks(ptr%ob_index)), &
!                        ' platform: ',trim(this%platforms(ptr%ob_index)), &
!                        ' lat: ',this%lats(ptr%ob_index), &
!                        ' lon: ',this%lons(ptr%ob_index), &
!                        ' elev: ',this%elevs(ptr%ob_index), &
!                        ' obs: ',this%obs(ptr%ob_index), &
!                        ' back: ',this%backs(ptr%ob_index)
!                   write(LDT_logunit,*) &
!                        '-----------------------------------------------------'
               else
                  reject_all = .true.
                  exit ! out of i loop
               end if
               ptr => ptr%next
            end do ! i

            if (reject_all) then
               ! We will reject all reports from this platform.
               ! Start with the "j" report, which is not stored in the linked
               ! list.
               this%qcs(j) = QC_REJECT
               this%reject_reasons(j) = "DUP_QC REJECT"
               total_reject_count = total_reject_count + 1
               ! write(LDT_logunit,*) &
!                     '[INFO] dupQC rejecting ob1 j: ', &
!                     j, &
!                     ' network: ',trim(this%networks(j)), &
!                     ' platform: ',trim(this%platforms(j)), &
!                     ' lat: ',this%lats(j), &
!                     ' lon: ',this%lons(j), &
!                     ' elev: ',this%elevs(j), &
!                     ' obs: ',this%obs(j), &
!                     ' back: ',this%backs(j)

               ! Now reject the duplicate reports in the linked list.
               ptr => head
               do i = 1, count_dups
                  this%qcs(ptr%ob_index) = QC_REJECT
                  this%reject_reasons(j) = "DUP_QC REJECT"
                  total_reject_count = total_reject_count + 1
!                   write(LDT_logunit,*) &
!                        '[INFO] dupQC rejecting ob1 i: ', &
!                        ptr%ob_index, &
!                        ' network: ',trim(this%networks(ptr%ob_index)), &
!                        ' platform: ',trim(this%platforms(ptr%ob_index)), &
!                        ' lat: ',this%lats(ptr%ob_index), &
!                        ' lon: ',this%lons(ptr%ob_index), &
!                        ' elev: ',this%elevs(ptr%ob_index), &
!                        ' obs: ',this%obs(ptr%ob_index), &
!                        ' back: ',this%backs(ptr%ob_index)

                  ptr => ptr%next
               end do ! i
               !write(LDT_logunit,*) &
               !     '------------------------------------------------------'
            end if ! reject_all
         end if ! count_dups .gt. 1 .and. .not. location_issue

         ! If we have exactly one duplicate: reject duplicate
         ! if it is an exact copy; otherwise, attempt superob.
         if (count_dups .eq. 1) then
            ptr => head
            diff = this%obs(ptr%ob_index) - this%obs(j)
            if (diff .eq. 0) then
               this%qcs(ptr%ob_index) = QC_REJECT
               this%reject_reasons(j) = "DUP_QC REJECT"
               total_reject_count = total_reject_count + 1
!                write(LDT_logunit,*) &
!                     '[INFO] dupQC rejecting exact duplicate2 ob j: ', &
!                     ptr%ob_index, &
!                     ' network: ',trim(this%networks(ptr%ob_index)), &
!                     ' platform: ',trim(this%platforms(ptr%ob_index)), &
!                     ' lat: ',this%lats(ptr%ob_index), &
!                     ' lon: ',this%lons(ptr%ob_index), &
!                     ' elev: ',this%elevs(ptr%ob_index), &
!                     ' obs: ',this%obs(ptr%ob_index), &
!                     ' back: ',this%backs(ptr%ob_index)
!                write(LDT_logunit,*) &
!                     '------------------------------------------------------'

            else if (diff*diff .gt. this%ob_err_vars(j)) then
               ! Obs are too different.  Reject both of them.
               this%qcs(j) = QC_REJECT
               this%reject_reasons(j) = "DUP_QC REJECT"
               total_reject_count = total_reject_count + 1
!                write(LDT_logunit,*) &
!                     '[INFO] dupQC rejecting2 ob j: ', &
!                     j, &
!                     ' network: ',trim(this%networks(j)), &
!                     ' platform: ',trim(this%platforms(j)), &
!                     ' lat: ',this%lats(j), &
!                     ' lon: ',this%lons(j), &
!                     ' elev: ',this%elevs(j), &
!                     ' obs: ',this%obs(j), &
!                     ' back: ',this%backs(j)

               this%qcs(ptr%ob_index) = QC_REJECT
               this%reject_reasons(j) = "DUP_QC REJECT"
               total_reject_count = total_reject_count + 1

!                write(LDT_logunit,*) &
!                     '[INFO] dupQC rejecting2 ob j: ', &
!                     ptr%ob_index, &
!                     ' network: ',trim(this%networks(ptr%ob_index)), &
!                     ' platform: ',trim(this%platforms(ptr%ob_index)), &
!                     ' lat: ',this%lats(ptr%ob_index), &
!                     ' lon: ',this%lons(ptr%ob_index), &
!                     ' elev: ',this%elevs(ptr%ob_index), &
!                     ' obs: ',this%obs(ptr%ob_index), &
!                     ' back: ',this%backs(ptr%ob_index)

               write(LDT_logunit,*) &
                    '------------------------------------------------------'

            else
               ! Create superob.
               mean = 0.5 * (this%obs(ptr%ob_index) + this%obs(j))

!                write(LDT_logunit,*) &
!                     '[INFO] dupQC will create superob from j: ', &
!                     j, &
!                     ' network: ',trim(this%networks(j)), &
!                     ' platform: ',trim(this%platforms(j)), &
!                     ' lat: ',this%lats(j), &
!                     ' lon: ',this%lons(j), &
!                     ' elev: ',this%elevs(j), &
!                     ' obs: ',this%obs(j), &
!                     ' back: ',this%backs(j)

!                write(LDT_logunit,*) &
!                     '[INFO] dupQC will create superob from i: ', &
!                     ptr%ob_index, &
!                     ' network: ',trim(this%networks(ptr%ob_index)), &
!                     ' platform: ',trim(this%platforms(ptr%ob_index)), &
!                     ' lat: ',this%lats(ptr%ob_index), &
!                     ' lon: ',this%lons(ptr%ob_index), &
!                     ' elev: ',this%elevs(ptr%ob_index), &
!                     ' obs: ',this%obs(ptr%ob_index), &
!                     ' back: ',this%backs(ptr%ob_index)

               network = trim(this%networks(j))
               platform = trim(this%platforms(j))
               back = this%backs(j)
               newlat = this%lats(j)
               newlon = this%lons(j)
               newelev = this%elevs(j)
               ob_err_var = this%ob_err_vars(j)

               ! write(LDT_logunit,*) &
!                     '[INFO] dupQC new superob is : ', &
!                     ' network: ',trim(this%networks(j)), &
!                     ' platform: ',trim(this%platforms(j)), &
!                     ' obs: ',mean
!                write(LDT_logunit,*) &
!                     '------------------------------------------------------'

               call this%append_ob(network, platform, mean, newlat, newlon,&
                    newelev, ob_err_var, back=back)

               total_create_count = total_create_count + 1

               ! Reject originals
               this%qcs(j) = QC_REJECT
               this%reject_reasons(j) = "DUP_QC MERGED"
               this%qcs(ptr%ob_index) = QC_REJECT
               this%reject_reasons(ptr%ob_index) = "DUP_QC MERGED"
               total_merge_count = total_merge_count + 2
            end if
         end if ! count_dups .eq. 1 .and. .not. location_issue

         ! Clean up linked list and move on
         do
            ptr => head
            if (associated(head%next)) then
               head => head%next
               deallocate(ptr)
            else
               nullify(head,tail)
               deallocate(ptr)
               nullify(ptr)
               exit ! Done with linked list
            end if
         end do

      end do ! j

      write(LDT_logunit,*) &
           '[INFO] dupQC rejected ',total_reject_count,' obs and merged ', &
           total_merge_count
      write(LDT_logunit,*) &
           '[INFO] dupQC created ',total_create_count,' super obs'

      if (total_reject_count + total_merge_count + total_create_count > 0) then
         call this%resort_bad_obs()
      end if

      call cpu_time(t2)
      write(LDT_logunit,*) &
           '[INFO] Elapsed time in dupQC is ',t2-t1,' seconds'

   end subroutine LDT_bratseth_run_dup_qc

   ! Reject observations over water
   subroutine LDT_bratseth_run_water_qc(this,n,ncols,nrows,landmask, &
        silent_rejects)

      ! Imports
      use LDT_coreMod, only: LDT_rc, LDT_domain
      use LDT_logMod, only: LDT_logunit
      use map_utils, only: latlon_to_ij

      ! Defaults
      implicit none

      ! Arguments
      class(LDT_bratseth_t), intent(inout) :: this
      integer, intent(in) :: n
      integer, intent(in) :: ncols
      integer, intent(in) :: nrows
      real, intent(in) :: landmask(ncols,nrows)
      logical,optional,intent(in) :: silent_rejects

      ! Local variables
      integer :: num_good_obs
      logical :: silent_rejects_local
      real :: t1, t2
      integer :: reject_count
      integer :: j
      real :: xpt, ypt
      integer :: c,r

      ! Sanity check
      num_good_obs = this%count_good_obs()
      if (num_good_obs .eq. 0) then
         write(LDT_logunit,*)&
              '[INFO] waterQC found no good observations to test'
         return
      end if

      silent_rejects_local = .false.
      if (present(silent_rejects)) then
         silent_rejects_local = silent_rejects
      end if

      call cpu_time(t1)

      reject_count = 0
      do j = 1, num_good_obs

         ! Find the column,row of this ob in the LDT domain
         call latlon_to_ij(LDT_domain(n)%ldtproj, &
              this%lats(j),this%lons(j), &
              xpt,ypt)
         c = nint(xpt)
         if (c > ncols) then
            c = c - ncols
         else if (c < 1) then
            c = c + ncols
         end if
         r = min(nrows,max(1,nint(ypt)))
         if (landmask(c,r) .lt. 0.5) then
            this%qcs(j) = QC_REJECT
            this%reject_reasons(j) = "WATER_QC REJECT"
            reject_count = reject_count + 1
            !if (.not. silent_rejects_local) then
            !   write(LDT_logunit,*) &
            !        '[INFO] waterQC rejecting observation i: ',j, &
            !        ' network: ',trim(this%networks(j)), &
            !        ' platform: ',trim(this%platforms(j)), &
            !        ' lat: ',this%lats(j), &
            !        ' lon: ',this%lons(j), &
            !        ' elev: ',this%elevs(j), &
            !        ' obs: ',this%obs(j)
            !end if
         end if
      end do ! j

      write(LDT_logunit,*)&
           '[INFO] waterQC rejected ',reject_count,' observations'

      if (reject_count > 0) then
         call this%resort_bad_obs()
      end if

      call cpu_time(t2)
      write(LDT_logunit,*) &
           '[INFO] Elapsed time in waterQC is ',t2-t1,' seconds'

   end subroutine LDT_bratseth_run_water_qc

   ! Creates "superobs" out of close observations.  Each close observation is
   ! first checked for unacceptable deviation from the mean of the close
   ! obs, and rejected if deviation is too high.  Superobs are considered
   ! "close" if they are in the same LIS grid box.  Based on Lespinas et al
   ! (2015).
   subroutine LDT_bratseth_run_superstat_qc(this, n, new_name, &
        ncols, nrows, silent_rejects)

      ! Imports
      use LDT_coreMod, only: LDT_domain, LDT_rc
      use LDT_logMod, only: LDT_logunit
      use map_utils, only: latlon_to_ij

      ! Defaults
      implicit none

      ! Arguments
      class(LDT_bratseth_t), intent(inout) :: this
      integer, intent(in) :: n
      character(len=32), intent(in) :: new_name
      integer, intent(in) :: ncols
      integer, intent(in) :: nrows
      logical, optional, intent(in) :: silent_rejects

      ! Local variables
      integer :: num_good_obs
      integer :: num_rejected_obs, num_merged_obs, num_superobs
      logical :: silent_rejects_local
      character(len=10) :: network_new
      character(len=32) :: platform_new
      integer, allocatable :: actions(:)
      real, allocatable :: means(:,:)
      real, allocatable :: superobs(:,:), superbacks(:,:)
      real, allocatable :: superlats(:,:),superlons(:,:),superelevs(:,:)
      real, allocatable :: superob_err_vars(:,:)
      integer, allocatable :: super_counts(:,:)
      integer :: icount
      real :: threshold
      integer :: c,r,j
      real :: xpt,ypt
      real :: t1, t2

      ! Sanity check
      num_good_obs = this%count_good_obs()
      if (num_good_obs .eq. 0) then
         write(LDT_logunit,*) &
              '[INFO] superstatQC found no good observations to test'
         return
      end if

      silent_rejects_local = .false.
      if (present(silent_rejects)) then
         silent_rejects_local = silent_rejects
      end if

      ! Determine which type of observations will be "superobbed"
      network_new = trim(new_name)
      platform_new = trim(new_name)

      call cpu_time(t1)

      allocate(actions(num_good_obs))
      actions(:) = 0

      allocate(means(ncols,nrows))
      means(:,:) = 0
      allocate(superobs(ncols,nrows))
      superobs(:,:) = 0
      allocate(superbacks(ncols,nrows))
      superbacks(:,:) = 0
      allocate(superlats(ncols,nrows))
      superlats(:,:) = 0
      allocate(superlons(ncols,nrows))
      superlons(:,:) = 0
      allocate(superelevs(ncols,nrows))
      superelevs(:,:) = 0
      allocate(superob_err_vars(ncols,nrows))
      superob_err_vars(:,:) = 0
      allocate(super_counts(ncols,nrows))
      super_counts(:,:) = 0

      ! Find all acceptable obs in each LDT grid box, and calculate mean
      ! values per box.
      do j = 1,num_good_obs

         ! Find the column,row of this ob in the LDT domain
         call latlon_to_ij(LDT_domain(n)%ldtproj, &
              this%lats(j),this%lons(j), &
              xpt,ypt)
         c = nint(xpt)
         if (c > ncols) then
            c = c - ncols
         else if (c < 1) then
            c = c + ncols
         end if
         r = min(nrows,max(1,nint(ypt)))

         ! Add contribution to local mean.
         means(c,r) = means(c,r) + this%obs(j)
         super_counts(c,r) = super_counts(c,r) + 1
      end do ! j

      do r = 1, LDT_rc%lnr(n)
         do c = 1, LDT_rc%lnc(n)
            if (super_counts(c,r) .gt. 1) then
               means(c,r) = means(c,r) / real(super_counts(c,r))
            end if
         end do ! c
      end do ! r

      ! The QC step:  Reject observations that deviate too much from local
      ! mean.
      do j = 1,num_good_obs

         ! Find the column,row of this ob in the LDT domain
         call latlon_to_ij(LDT_domain(n)%ldtproj, &
              this%lats(j),this%lons(j), &
              xpt,ypt)
         c = nint(xpt)
         if (c > ncols) then
            c = c - ncols
         else if (c < 1) then
            c = c + ncols
         end if
         r = min(nrows,max(1,nint(ypt)))

         ! Need at least two local obs for this test to make sense.
         if (super_counts(c,r) .lt. 2) cycle

         icount = super_counts(c,r)
         threshold = 3 * this%ob_err_vars(j) * &
              sqrt(real(icount) / real(icount - 1))

         if (abs(means(c,r) - this%obs(j)) .gt. threshold) then
            actions(j) = -1 ! Reject
         else
            actions(j) =  1 ! Consider for superob
         end if

      end do ! j

      ! Calculate superobs from remaining good observations.
      super_counts(:,:) = 0 ! Reset since some obs were rejected above.
      do j = 1,num_good_obs

         ! Only use the obs that passed the test above.
         if (actions(j) .ne. 1) cycle

         ! Find the column,row of this ob in the LDT domain
         call latlon_to_ij(LDT_domain(n)%ldtproj, &
              this%lats(j),this%lons(j), &
              xpt,ypt)
         c = nint(xpt)
         if (c > ncols) then
            c = c - ncols
         else if (c < 1) then
            c = c + ncols
         end if
         r = min(nrows,max(1,nint(ypt)))

         super_counts(c,r) = super_counts(c,r) + 1
         superobs(c,r) = superobs(c,r) + this%obs(j)
         superbacks(c,r) = superbacks(c,r) + this%backs(j)
         superlats(c,r) = superlats(c,r) + this%lats(j)
         superlons(c,r) = superlons(c,r) + this%lons(j)
         superelevs(c,r) = superelevs(c,r) + this%elevs(j)
         superob_err_vars(c,r) = superob_err_vars(c,r) + this%ob_err_vars(j)

      end do ! j

      do r = 1, LDT_rc%lnr(n)
         do c = 1, LDT_rc%lnc(n)

            ! Need at least two obs to merge
            if (super_counts(c,r) .lt. 2) cycle

            superobs(c,r) = superobs(c,r) / real(super_counts(c,r))
            superbacks(c,r) = superbacks(c,r) / real(super_counts(c,r))
            superlats(c,r) = superlats(c,r) / real(super_counts(c,r))
            superlons(c,r) = superlons(c,r) / real(super_counts(c,r))
            superelevs(c,r) = superelevs(c,r) / real(super_counts(c,r))
            superob_err_vars(c,r) = superob_err_vars(c,r) / &
                 real(super_counts(c,r))

         end do ! c
      end do ! r

      ! Make sure to preserve obs that passed the QC test above but were not
      ! merged (not enough good obs in the grid box).
      do j = 1,num_good_obs

         ! Only use the obs that passed the test above.
         if (actions(j) .ne. 1) cycle

         ! Find the column,row of this ob in the LDT domain
         call latlon_to_ij(LDT_domain(n)%ldtproj, &
              this%lats(j),this%lons(j), &
              xpt,ypt)
         c = nint(xpt)
         if (c > ncols) then
            c = c - ncols
         else if (c < 1) then
            c = c + ncols
         end if
         r = min(nrows,max(1,nint(ypt)))

         if (super_counts(c,r) .lt. 2) then
            actions(j) = 0 ! Preserve this ob
         end if
      end do ! j

      ! Now update the observation object with the rejected obs and merged
      ! obs
      num_merged_obs = 0
      num_rejected_obs = 0
      do j = 1, num_good_obs
         if (actions(j) .eq. -1) then
            this%qcs(j) = QC_REJECT
            this%reject_reasons(j) = "SUPERSTAT_QC REJECT"
            ! if (.not. silent_rejects_local) then
!                write(LDT_logunit,*) &
!                     '[INFO] superstatQC rejection j: ',j, &
!                     ' net: ',trim(this%networks(j)), &
!                     ' platform: ',trim(this%platforms(j)), &
!                     ' lat: ',this%lats(j), &
!                     ' lon: ',this%lons(j), &
!                     ' elev: ',this%elevs(j), &
!                     ' obs: ',this%obs(j), &
!                     ' back: ',this%backs(j)
!             end if
            num_rejected_obs = num_rejected_obs + 1
         else if (actions(j) .eq. 1) then
            this%qcs(j) = QC_REJECT ! Was merged into superob
            this%reject_reasons(j) = "SUPERSTAT_QC MERGED"
!             if (.not. silent_rejects_local) then
!                write(LDT_logunit,*) &
!                     '[INFO] superstatQC will create superob from j: ',j, &
!                     ' net: ',trim(this%networks(j)), &
!                     ' platform: ',trim(this%platforms(j)), &
!                     ' lat: ',this%lats(j), &
!                     ' lon: ',this%lons(j), &
!                     ' elev: ',this%elevs(j), &
!                     ' obs: ',this%obs(j), &
!                     ' back: ',this%backs(j)
!             end if
            num_merged_obs = num_merged_obs + 1
         end if
      end do ! j

      write(LDT_logunit,*) &
           '[INFO] superstatQC rejected ',num_rejected_obs, &
           ' obs and merged ',num_merged_obs,' obs'

      ! Finally, add the superobs to the object
      num_superobs = 0
      do r = 1, LDT_rc%lnr(n)
         do c = 1, LDT_rc%lnc(n)
            if (super_counts(c,r) .lt. 2) cycle

            call this%append_ob(network_new, platform_new, &
                 superobs(c,r), &
                 superlats(c,r), &
                 superlons(c,r), &
                 superelevs(c,r), &
                 superob_err_vars(c,r), &
                 superbacks(c,r))

            num_superobs = num_superobs + 1
         end do ! c
      end do ! r

      write(LDT_logunit,*) &
           '[INFO] superstatQC created ',num_superobs,' super obs'

      ! Clean up
      deallocate(means)
      deallocate(actions)
      deallocate(superobs)
      deallocate(superbacks)
      deallocate(superlats)
      deallocate(superlons)
      deallocate(superelevs)
      deallocate(superob_err_vars)
      deallocate(super_counts)

      if (num_rejected_obs + num_merged_obs + num_superobs > 0) then
         call this%resort_bad_obs()
      end if

      call cpu_time(t2)
      write(LDT_logunit,*) &
           '[INFO] Elapsed time in superstatQC is ',t2-t1,' seconds'

   end subroutine LDT_bratseth_run_superstat_qc

   ! Reject obs that differ "too much" from background field.  Threshold based
   ! on sum of observation and background error variances.  From Lopez (2013).
   ! This assumes both background and observations are unbiased, and large
   ! difference implies gross error in observation.
   subroutine LDT_bratseth_run_back_qc(this,back_err_var,silent_rejects)

      ! Imports
      use LDT_logMod, only: LDT_logunit

      ! Defaults
      implicit none

      ! Arguments
      class(LDT_bratseth_t), intent(inout) :: this
      real, intent(in) :: back_err_var
      logical, optional, intent(in) :: silent_rejects

       ! Local variables
      integer :: num_good_obs
      logical :: silent_rejects_local
      real :: t1, t2
      integer :: reject_count
      integer :: j
      real :: threshold, abs_diff

      ! Sanity check
      num_good_obs = this%count_good_obs()
      if (num_good_obs .eq. 0) then
         write(LDT_logunit,*)&
              '[INFO] backQC found no good observations to test'
         return
      endif

      silent_rejects_local = .false.
      if (present(silent_rejects)) then
         silent_rejects_local = silent_rejects
      end if

      call cpu_time(t1)

      reject_count = 0
      do j = 1, num_good_obs

         threshold = 5*sqrt(back_err_var + this%ob_err_vars(j))
         abs_diff = abs(this%obs(j) - this%backs(j))

         if (abs_diff .gt. threshold) then
            this%qcs(j) = QC_REJECT
            this%reject_reasons(j) = "BACK_QC REJECT"

            reject_count = reject_count + 1

!             if (.not. silent_rejects_local) then
!                write(LDT_logunit,*) &
!                     '[INFO] backQC rejecting observation j: ',j, &
!                     ' network: ',trim(this%networks(j)), &
!                     ' platform: ',trim(this%platforms(j)), &
!                     ' lat: ',this%lats(j), &
!                     ' lon: ',this%lons(j), &
!                     ' elev: ',this%elevs(j), &
!                     ' obs: ',this%obs(j), &
!                     ' back: ',this%backs(j), &
!                     ' abs diff: ', abs(this%obs(j) - this%backs(j)), &
!                     ' threshold ', threshold
!             end if
         end if

      end do ! j

      write(LDT_logunit,*)&
           '[INFO] backQC rejected ',reject_count,' observations'

      if (reject_count > 0) then
         call this%resort_bad_obs()
      end if

      call cpu_time(t2)
      write(LDT_logunit,*) &
           '[INFO] Elapsed time in backQC is ',t2-t1,' seconds'

   end subroutine LDT_bratseth_run_back_qc

   ! Skewed error check against the background.  If the observation indicates
   ! snow depth at least X less than the background, reject the observation.
   ! From Brasnett (1999), where X is usually 40 cm.
   subroutine LDT_bratseth_run_skewed_back_qc(this,threshold)

      ! Imports
      use LDT_logMod, only: LDT_logunit

      ! Defaults
      implicit none

      ! Arguments
      class(LDT_bratseth_t), intent(inout) :: this
      real, intent(in) :: threshold

      ! Local variables
      integer :: num_good_obs
      real :: diff
      integer :: j
      real :: t1, t2
      integer :: total_reject_count

      ! Sanity check
      num_good_obs = this%count_good_obs()
      if (num_good_obs .eq. 0) then
         write(LDT_logunit,*) &
              '[INFO] skewedBackQC found no good observations to test'
         return
      end if

      call cpu_time(t1)

      total_reject_count = 0

      ! Compare each good ob to the background, rejecting those that are
      ! grossly too small
      do j = 1, num_good_obs

         diff = this%backs(j) - this%obs(j)
         if (diff .gt. threshold) then
            this%qcs(j) = QC_REJECT
            this%reject_reasons(j) = "SKEWED_BACK_QC REJECT"
            total_reject_count = total_reject_count + 1
            !write(LDT_logunit,*) &
            !     '[INFO] skewedBackQC rejecting ob j: ', j, &
            !     ' net: ',trim(this%networks(j)), &
            !     ' platform: ',trim(this%platforms(j)), &
            !     ' lat: ',this%lats(j), &
            !     ' lon: ',this%lons(j), &
            !     ' elev: ',this%elevs(j), &
            !     ' obs: ',this%obs(j), &
            !     ' back: ',this%backs(j)
            !write(LDT_logunit,*) &
            !           '------------------------------------------------------'
         end if
      end do ! j

      write(LDT_logunit,*) &
           '[INFO] skewedBackQC rejected ',total_reject_count,' obs'

      if (total_reject_count > 0) then
         call this%resort_bad_obs()
      end if

      call cpu_time(t2)
      write(LDT_logunit,*) &
           '[INFO] Elapsed time in skewedBackQC is ',t2-t1,' seconds'

   end subroutine LDT_bratseth_run_skewed_back_qc

   ! Compare reported station elevation to interpolated LDT value, and
   ! reject ob if the elevations differ too much.  Based on Brasnett (1999).
   subroutine LDT_bratseth_run_elev_qc(this,n,ncols,nrows,elevations,threshold)

      ! Imports
      use LDT_coreMod, only: LDT_domain
      use LDT_logMod, only: LDT_logunit
      use map_utils

      ! Defaults
      implicit none

      ! Arguments
      class(LDT_bratseth_t), intent(inout) :: this
      integer,intent(in) :: n
      integer,intent(in) :: ncols
      integer,intent(in) :: nrows
      real, intent(in) :: elevations(ncols,nrows)
      real, intent(in) :: threshold

      ! Local variables
      integer :: num_good_obs
      real :: absdiff
      integer :: c,r,j
      real :: xpt,ypt
      real :: t1, t2
      integer :: total_reject_count

      ! Sanity check
      num_good_obs = this%count_good_obs()
      if (num_good_obs .eq. 0) then
         write(LDT_logunit,*) &
              '[INFO] elevQC found no good observations to test'
         return
      end if

      call cpu_time(t1)

      total_reject_count = 0

      ! Compare elevation of each good ob to the interpolated terrain value,
      ! rejecting those that differ too much.
      do j = 1,num_good_obs

         ! Find the LDT grid box containing the observation.
         call latlon_to_ij(LDT_domain(n)%ldtproj, &
              this%lats(j),this%lons(j), &
              xpt,ypt)
         c = nint(xpt)
         if (c > ncols) then
            c = c - ncols
         else if (c < 1) then
            c = c + ncols
         end if
         r = min(nrows,max(1,nint(ypt)))

         absdiff = abs(elevations(c,r) - this%elevs(j))
         if (absdiff .gt. threshold) then
            this%qcs(j) = QC_REJECT
            this%reject_reasons(j) = "ELEV_QC REJECT"
            total_reject_count = total_reject_count + 1
            !write(LDT_logunit,*) &
            !     '[INFO] elevQC rejecting ob j: ', j, &
            !     ' net: ',trim(this%networks(j)), &
            !     ' platform: ',trim(this%platforms(j)), &
            !     ' lat: ',this%lats(j), &
            !     ' lon: ',this%lons(j), &
            !     ' elev: ',this%elevs(j), &
            !     ' LDT terrain: ',elevations(c,r), &
            !     ' obs: ',this%obs(j), &
            !     ' back: ',this%backs(j)
            !write(LDT_logunit,*) &
            !           '------------------------------------------------------'
         end if
      end do ! j

      write(LDT_logunit,*) &
           '[INFO] elevQC rejected ',total_reject_count,' obs'

      if (total_reject_count > 0) then
         call this%resort_bad_obs()
      end if

      call cpu_time(t2)
      write(LDT_logunit,*) &
           '[INFO] Elapsed time in elevQC is ',t2-t1,' seconds'

   end subroutine LDT_bratseth_run_elev_qc

   ! Get lat/lon of observation
   subroutine LDT_bratseth_get_lat_lon(this,job,lat,lon)
      implicit none
      class(LDT_bratseth_t), intent(in) :: this
      integer, intent(in) :: job
      real, intent(out) :: lat
      real, intent(out) :: lon
      lat = this%lats(job)
      lon = this%lons(job)
   end subroutine LDT_bratseth_get_lat_lon

   ! Set the interpolated background value for a specific ob
   subroutine LDT_bratseth_set_back(this,job,back)
      implicit none
      class(LDT_bratseth_t), intent(inout) :: this
      integer, intent(in) :: job
      real, intent(in) :: back
      this%backs(job) = back
   end subroutine LDT_bratseth_set_back

   ! Reject specified observation.  Useful if background is missing or
   ! snomask disagrees.
   subroutine LDT_bratseth_reject_ob(this,job,reason)
      implicit none
      class(LDT_bratseth_t), intent(inout) :: this
      integer, intent(in) :: job
      character(len=*), intent(in) :: reason
      this%qcs(job) = QC_REJECT
      this%reject_reasons(job) = trim(reason)
   end subroutine LDT_bratseth_reject_ob

   ! Moves rejected obs towards end of Bratseth arrays
   subroutine LDT_bratseth_resort_bad_obs(this)

      ! Imports
      use LDT_logMod, only: LDT_logunit, LDT_endrun

      ! Defaults
      implicit none

      ! Arguments
      class(LDT_bratseth_t), intent(inout) :: this

      ! Local variables
      integer :: nobs
      integer :: iob,job
      logical :: found_replacement
      character*10 :: network
      character*32 :: platform
      real :: ob, lat, lon, elev, ob_err_var, back
      integer :: qc
      character*80 :: failed_qc_test

      ! Sanity check
      nobs = this%nobs
      if (nobs .eq. 0) return

      do job = 1, nobs
         if (this%qcs(job) == QC_REJECT) then

            ! Look for a good ob further down the list
            found_replacement = .false.
            do iob = job+1,nobs
               if (this%qcs(iob) .ne. QC_REJECT) then

                  ! Swap the good and bad obs
                  found_replacement = .true.

                  network = this%networks(job)
                  platform = this%platforms(job)
                  ob = this%obs(job)
                  lat = this%lats(job)
                  lon = this%lons(job)
                  elev = this%elevs(job)
                  ob_err_var = this%ob_err_vars(job)
                  back = this%backs(job)
                  qc = this%qcs(job)
                  failed_qc_test = this%reject_reasons(job)

                  this%networks(job) = this%networks(iob)
                  this%platforms(job) = this%platforms(iob)
                  this%obs(job) = this%obs(iob)
                  this%lats(job) = this%lats(iob)
                  this%lons(job) = this%lons(iob)
                  this%elevs(job) = this%elevs(iob)
                  this%ob_err_vars(job) = this%ob_err_vars(iob)
                  this%backs(job) = this%backs(iob)
                  this%qcs(job) = this%qcs(iob)
                  this%reject_reasons(job) = this%reject_reasons(iob)

                  this%networks(iob) = network
                  this%platforms(iob) = platform
                  this%obs(iob) = ob
                  this%lats(iob) = lat
                  this%lons(iob) = lon
                  this%elevs(iob) = elev
                  this%ob_err_vars(iob) = ob_err_var
                  this%backs(iob) = back
                  this%qcs(iob) = qc
                  this%reject_reasons(iob) = failed_qc_test

                  exit ! Get out of iob loop
               end if
            end do ! iob

            ! If we didn't find a replacement, all the remaining obs are bad.
            if (.not. found_replacement) then
               exit
            end if

         end if
      end do ! job

   end subroutine LDT_bratseth_resort_bad_obs

   ! Return observed value
   subroutine LDT_bratseth_get_ob_value(this,job,ob_value)
      implicit none
      class(LDT_bratseth_t), intent(in) :: this
      integer, intent(in) :: job
      real, intent(out) :: ob_value
      ob_value = this%obs(job)
   end subroutine LDT_bratseth_get_ob_value

   ! Reject obs where snow is not permitted
   subroutine LDT_bratseth_run_nosnow_qc(this,nc,nr,snow_poss)

      ! Imports
      use LDT_coreMod, only: LDT_domain, LDT_rc
      use LDT_logMod, only: LDT_logunit
      use map_utils

      ! Defaults
      implicit none

      ! Arguments
      class(LDT_bratseth_t), intent(inout) :: this
      integer, intent(in) :: nc
      integer, intent(in) :: nr
      integer*1, intent(in) :: snow_poss(nc,nr)

      ! Local variables
      integer :: num_good_obs
      real :: t1,t2
      real :: rc,rr
      integer :: c, r
      integer :: reject_count
      integer :: j

      ! Sanity check
      num_good_obs = this%count_good_obs()
      if (num_good_obs .eq. 0) then
         write(LDT_logunit,*)&
              '[INFO] nosnowQC found good observations to test'
         return
      endif

      call cpu_time(t1)

      reject_count = 0

      do j = 1, num_good_obs

         ! EMK Only reject the ob when snow is not possible in low- and &
         ! mid-latitudes.  At high latitudes, we assume snow is possible.
         if (this%lats(j) <= -40.0) cycle
         if (this%lats(j) >=  40.0) cycle

         call latlon_to_ij(LDT_domain(1)%ldtproj, &
              this%lats(j), &
              this%lons(j), &
              rc, &
              rr)

         c = nint(rc)
         if (c < 1) then
            c = c + LDT_rc%lnc(1)
         else if (c > nc) then
            c = c - LDT_rc%lnc(1)
         end if
         r = min(LDT_rc%lnr(1),max(1,nint(rr)))
         if (snow_poss(c,r) == 0) then
            !write(LDT_logunit,*) &
            !     '[INFO] REJECTING OB FOR BEING IN NO-SNOW REGION: ', &
            !     ' network: ',trim(this%networks(j)), &
            !     ' platform: ',trim(this%platforms(j)), &
            !     ' lat: ',this%lats(j), &
            !     ' lon: ',this%lons(j), &
            !     ' elev: ',this%elevs(j), &
            !     ' obs: ',this%obs(j)
            this%qcs(j) = QC_REJECT
            this%reject_reasons(j) = "NOSNOW_QC REJECT"
            reject_count = &
                 reject_count + 1
         end if
      end do ! j

      write(LDT_logunit,*)&
           '[INFO] nosnowQC rejected ',reject_count,' observations'

      if (reject_count > 0) then
         call this%resort_bad_obs()
      end if

      call cpu_time(t2)
      write(LDT_logunit,*) &
           '[INFO] Elapsed time in nosnowQC is ',t2-t1,' seconds'

   end subroutine LDT_bratseth_run_nosnow_qc

   ! Reject obs that are missing elevations
   subroutine LDT_bratseth_run_missing_elev_qc(this)

      ! Imports
      use LDT_logMod, only: LDT_logunit

      ! Defaults
      implicit none

      ! Arguments
      class(LDT_bratseth_t), intent(inout) :: this

      ! Local variables
      real :: t1, t2
      integer :: reject_count, num_good_obs
      integer :: j

      ! Sanity check
      num_good_obs = this%count_good_obs()
      if (num_good_obs .eq. 0) then
         write(LDT_logunit,*)&
              '[INFO] missingElevQC found no good observations to test'
         return
      endif

      call cpu_time(t1)

      reject_count = 0

      do j = 1, num_good_obs
         if (this%elevs(j) == -1000) then
            !write(LDT_logunit,*) &
            !     '[INFO] REJECTED OB FOR MISSING ELEVATION: ', &
            !     ' network: ',trim(this%networks(j)), &
            !     ' platform: ',trim(this%platforms(j)), &
            !     ' lat: ',this%lats(j), &
            !     ' lon: ',this%lons(j), &
            !     ' elev: ',this%elevs(j), &
            !     ' obs: ',this%obs(j)
            this%qcs(j) = QC_REJECT
            this%reject_reasons(j) = 'MISSING_ELEV_QC REJECT'
            reject_count = &
                 reject_count + 1
         end if
      end do ! j

      write(LDT_logunit,*)&
           '[INFO] missingElevQC rejected ',reject_count,' observations'

      if (reject_count > 0) then
         call this%resort_bad_obs()
      end if

      call cpu_time(t2)
      write(LDT_logunit,*) &
           '[INFO] Elapsed time in missingElevQC is ',t2-t1,' seconds'

   end subroutine LDT_bratseth_run_missing_elev_qc

   ! Sort observations by ID.  Based on legacy SNODEP subroutine COMBSORT
   ! See https://en.wikipedia.org/wiki/Comb_sort
   subroutine LDT_bratseth_sort_obs_by_id(this)

      ! Imports
      use LDT_logmod, only : LDT_logunit

      ! Defaults
      implicit none

      ! Arguments
      class(LDT_bratseth_t), intent(inout) :: this

      ! Local variables
      integer :: gap
      integer :: switch
      character*10 :: ctemp10
      character*32 :: ctemp32
      character*80 :: ctemp80
      real :: rtemp
      integer :: itemp
      integer :: i,j

      ! Sanity check
      if (this%nobs == 0) return

      ! Initializations
      gap = this%nobs

      do
         switch = 0
         gap = max(int(float(gap) / 1.3), 1)
         if ( (gap == 9) .or. (gap == 10)) then
            gap = 11
         end if

         do i = 1, this%nobs - gap
            j = i + gap
            if (this%platforms(i) > this%platforms(j)) then

               ctemp10 = this%networks(i)
               this%networks(i) = this%networks(j)
               this%networks(j) = ctemp10

               ctemp32 = this%platforms(i)
               this%platforms(i) = this%platforms(j)
               this%platforms(j) = ctemp32

               rtemp = this%obs(i)
               this%obs(i) = this%obs(j)
               this%obs(j) = rtemp

               rtemp = this%lats(i)
               this%lats(i) = this%lats(j)
               this%lats(j) = rtemp

               rtemp = this%lons(i)
               this%lons(i) = this%lons(j)
               this%lons(j) = rtemp

               rtemp = this%elevs(i)
               this%elevs(i) = this%elevs(j)
               this%elevs(j) = rtemp

               rtemp = this%ob_err_vars(i)
               this%ob_err_vars(i) = this%ob_err_vars(j)
               this%ob_err_vars(j) = rtemp

               rtemp = this%backs(i)
               this%backs(i) = this%backs(j)
               this%backs(j) = rtemp

               itemp = this%qcs(i)
               this%qcs(i) = this%qcs(j)
               this%qcs(j) = itemp

               ctemp80 = this%reject_reasons(i)
               this%reject_reasons(i) = this%reject_reasons(j)
               this%reject_reasons(j) = ctemp80

               rtemp = this%inv_data_dens(i)
               this%inv_data_dens(i) = this%inv_data_dens(j)
               this%inv_data_dens(j) = rtemp

               rtemp = this%sum_ob_ests(i)
               this%sum_ob_ests(i) = this%sum_ob_ests(j)
               this%sum_ob_ests(j) = rtemp

               rtemp = this%anas(i)
               this%anas(i) = this%anas(j)
               this%anas(j) = rtemp

               switch = switch + 1
            end if
         end do ! i

         if ((switch > 0) .or. (gap > 1)) cycle
         exit

      end do

   end subroutine LDT_bratseth_sort_obs_by_id

   ! Print out any station with a snow depth
   subroutine LDT_bratseth_print_snowdepths(this,minprt)
      use LDT_logMod, only: LDT_logunit
      implicit none
      class(LDT_bratseth_t), intent(in) :: this
      real, intent(in) :: minprt
      integer :: i
      if (this%nobs == 0) return
      do i = 1, this%nobs
         if (this%obs(i) >= minprt) then
            write(LDT_logunit,7000) this%networks(i), &
                 this%platforms(i), this%lats(i), this%lons(i), &
                 int(this%elevs(i)), this%obs(i), trim(this%reject_reasons(i))
         end if
      end do ! i
7000  format (/,'[INFO] NETID = ',A5,' STATION ID = ',A32, &
           '   LAT = ',F6.2,'   LON = ',F7.2, &
           '   ELEV(M) = ',I5,'   DEPTH(M) = ', F7.5, &
           '   QC VERDICT = ',A)
   end subroutine LDT_bratseth_print_snowdepths

   ! Retrieve total obs
   function LDT_bratseth_count_all_obs(this) result (count)
      implicit none
      class(LDT_bratseth_t), intent(in) :: this
      integer :: count
      count = this%nobs
   end function LDT_bratseth_count_all_obs

   ! Mark good obs
   subroutine LDT_bratseth_mark_good_obs(this)
      implicit none
      class(LDT_bratseth_t), intent(inout) :: this
      integer :: i, num_good_obs
      num_good_obs = this%count_good_obs()
      do i = 1, num_good_obs
         this%reject_reasons(i) = "GOOD"
      end do ! i
   end subroutine LDT_bratseth_mark_good_obs
end module LDT_bratsethMod

