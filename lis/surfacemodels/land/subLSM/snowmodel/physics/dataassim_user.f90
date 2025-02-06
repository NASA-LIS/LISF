! dataassim_user.f90

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! subroutine DATAASSIM_USER does the SWE assimilation.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine DATAASSIM_USER(nx,ny,icorr_factor_index,&
     &  corr_factor,max_iter,deltax,deltay,xmn,ymn,nobs_dates,&
     &  print_inc,iday_init,imonth_init,iyear_init,dt,&
     &  output_path_wo_assim,xhour_init)

! Perform the required correction (precipitation and melt) factor
!   calculations.  To run the data assimilation routines requires
!   additional model inputs to be added here (inputs in addition
!   to those in the snowmodel.par file).  Then you must recompile
!   the code before running it.

! This program works when you have data at one or many points
!   (many individual grid cells), for one or more times.  And for
!   data (e.g., average swe) over areas of grid cells; there can be
!   many of these areas, at many different times.

! The input (swe) data file is assumed to be in the directory:
!   /snowmodel/swe_assim/.  See below for the details of how this
!   file is to be configured.

! The data assimilation produces a couple of extra files the user
!   can look at to study what the model did with the input data
!   file (fname_sweobs) that was provided.  First, there is a text
!   file written to /snowmodel/swe_assim/ that provides a summary
!   of the calculations that were done (see unit=77 in the code
!   below for what is written to this file).  Second, there is a
!   copy of the precipitation correction surfaces that were
!   generated.  This is also placed in /snowmodel/swe_assim/, the
!   file is a GrADS file (called corr_factor.gdat), and there is
!   also a corr_factor.ctl there that can be modified to fit the
!   current simulation.  The data layers in corr_factor.gdat arec
!   structured as follows:
!   The number of layers equals the total number of observation
!   times that were assimilated, plus the number of  years in the
!   assimilation.  The order is: the correction surface (cf) for
!   the time between the start of the simulation and observation
!   time 1 in year 1, the cf for the time between obs time 1 in
!   year 1 and obs time 2 in year 1, etc., then a cf==1.0 for the
!   time between the last obs time and the end of year 1 (or the
!   end of the simulation for a 1-year run).  Then this order
!   repeats for the rest of the years of the simulation.  In the
!   GrADS control file (corr_factor.ctl) these layers are assumed
!   to correspond to different times in the data file (although
!   the actual time increment defined in the .ctl file is not
!   really relevant: for example, t=1 corresponds to the cf for
!   obs 1, t=2 is for obs 2, t=3 is for the time between the last
!   obs time and the end of year 1, t=4 is for the obs 1 in year
!   2, etc.).

      use snowmodel_inc
      implicit none

      real deltax,deltay,beta,areas_flag,print_inc,dt,xhour_init
      real sprec_ratio(max_obs_dates),smelt_ratio(max_obs_dates)
      real corr_factor(nx_max,ny_max,max_obs_dates)
      real areas_mask(nx_max,ny_max)

      double precision xmn,ymn
      double precision xstn(nx_max*ny_max),ystn(nx_max*ny_max)

      integer icorr_factor_index(max_time_steps)
      integer iobs_rec(max_obs_dates)
      integer nobs_dates,nx,ny,max_iter,local_assim_flag,iday_init,&
     &  imonth_init,iyear_init,nyear,nyears,nobs_total,nobs_total_cfi

      character*89 fname_swed,fname_sspr,fname_ssmt
      character*80 fname_sweobs
      character*80 fname_sweobs_barnes_mask

      character*80 output_path_wo_assim
      integer trailing_blanks,i_len_wo
      integer i,j

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! BEGIN USER EDIT SECTION.
! BEGIN USER EDIT SECTION.

! Define how many years there are in this simulation.  Note that
!   these are SnowModel simulation years.  So, for example, year
!   1 might be 1 Sep 2007 - 31 Aug 2008, and year 2 might be
!   1 Sep 2008 - 31 Aug 2009, etc.
! NOTE: this is now put in the first line of the swe assimilation
!   input file (in 'fname_sweobs' below).
!     nyears = 20
!     nyears = 1

! A single data file describes the observation information that
!   will be used in the data assimilation.  This file contains the
!   following information in the following format.  The id can be
!   any number, x and y are easting and northing in m, and swe is
!   in m).

!   total_number_of_observation_dates_for_this_year
!   iyr imo idy (for this observation date)
!   number_of_stations_for_this_observation_date
!   id x y swe
!   id x y swe
!   id x y swe
!   iyr imo idy (for this observation date)
!   number_of_stations_for_this_observation_date
!   id x y swe
!   id x y swe

! For example:

!   2
!   2014 3 15
!   3
!   101 3456.7 23677.4 0.42
!   102 3556.3 25079.3 0.52
!   103 3106.2 29089.3 0.59
!   2014 4 1
!   2
!   101 3456.7 23677.4 0.48
!   103 3106.2 29089.3 0.62

! Then this repeats for each year of the assimilation (the input
!   file looks like the above for a single-year run, and for
!   multi-year runs the data for each following year is just
!   stacked on top of the previous year).  The example below is
!   for a two-year assimilation run.
! NOTE: the first '2=nyears' is because this is two-year run (if
!   this were a 1-year run, this would be a 1, etc.).

!   2
!   2
!   2014 3 15
!   3
!   101 3456.7 23677.4 0.42
!   102 3556.3 25079.3 0.52
!   103 3106.2 29089.3 0.59
!   2014 4 1
!   2
!   101 3456.7 23677.4 0.48
!   103 3106.2 29089.3 0.62
!   1
!   2015 3 25
!   2
!   101 3456.7 23677.4 0.23
!   102 3556.3 25079.3 0.32

! For the run where you have years with no data to assimilate,
!   and some years with data to assimilate, the code still
!   requires a "total_number_of_observation_dates_for_this_year"
!   line for each year.  So, if you have a year with no data, this
!   must be set to be 0 (zero).  If this is 0, the code sets the
!   correction factor to equal 1.0 for that simulation year, and
!   no adjustments are made to the precipitation for that year.
!   For example, if the first two, and fourth, years of a five-
!   year run have no data to assimilate, your input file would
!   look like:

!   5
!   0
!   0
!   1
!   2014 3 15
!   3
!   101 3456.7 23677.4 0.42
!   102 3556.3 25079.3 0.52
!   103 3106.2 29089.3 0.59
!   0
!   1
!   2015 3 25
!   2
!   101 3456.7 23677.4 0.23
!   102 3556.3 25079.3 0.32

! Provide the name of the data file that contains the observed swe
!   information (as described above).  This can no longer be
!   changed (defining it this way allows the error checking to be
!   done).
      fname_sweobs = 'swe_assim/swe_obs.dat'

! Define the file names of the swe depth (swed), annual summed snow
!   precipitation (sspr), and annual summed snowmelt (ssmt) outputs
!   from the first iteration of the data assimilation run.  In this
!   implementation of the data assimilation code, I have assumed
!   that the output files are those created by outputs_user.f,
!   where there is an individual file for each variable.
! NOTE: in the latest code version these paths have already been
!   defined in the snowmodel.par file.
!     fname_swed = 'outputs/wo_assim/swed.gdat'
!     fname_sspr = 'outputs/wo_assim/sspr.gdat'
!     fname_ssmt = 'outputs/wo_assim/ssmt.gdat'
      i_len_wo = 80 - trailing_blanks(output_path_wo_assim)
      fname_swed = output_path_wo_assim(1:i_len_wo)//'swed.gdat'
      fname_sspr = output_path_wo_assim(1:i_len_wo)//'sspr.gdat'
      fname_ssmt = output_path_wo_assim(1:i_len_wo)//'ssmt.gdat'

! THE PARAMETERS BELOW ARE RARELY CHANGED, UNLESS YOU ARE DOING AN
!   AREAS ASSIMILATION (INSTEAD OF ASSIMILATING POINT DATA).

! Beta controls the interpolation distance weights.  Beta = 1.0
!   will give you a very smooth field, and correction factor
!   distributions that may not produce swe's that exactly match
!   the observations.  Beta << 1.0 will give you correction factor
!   fields that go right through the data.  If you just have one
!   data point/area, beta is not used.
      beta = 1.0
!     beta = 0.1
!     beta = 0.5

! Define whether this simulation will be processing areas (data
!   within groups of grid cells: areas_flag = 1.0), or points
!   (single grid cells: areas_flag = 0.0).  Note that if you have
!   a combination of areas and points, you have to use the areas
!   option and treat each point like a single-grid-cell (small)
!   area.
      areas_flag = 0.0
!     areas_flag = 1.0

! If this is an areas simulation, open and read in the areas mask
!   data.  Note that here I assume that the area mask is a nx by ny
!   file with undef values everywhere except at the area 'stations'.
!   And that each 'station' area is given a 1.0, 2.0, etc. that
!   corresponds to the order of the station listing in the 'station'
!   data input file (the first 'station' listed has mask value = 1.0,
!   the second listed has mask value = 2.0, etc.
      if (areas_flag.eq.1.0) then
!       open(63,file=
!    &    '../1_topo_vege/8_mk_glac_front_mask/seals_fjord_mask.gdat',
!    &    form='unformatted',access='direct',recl=4*nx*ny)
!       read(63,rec=1) ((areas_mask(i,j),i=1,nx),j=1,ny)
! If you have two masks for two different observation dates, then
!   do something like the following.
!       open(63,file='swe_assim/zack_obs_mask.gdat',
!    &    form='unformatted',access='direct',recl=4*nx*ny)
!       read(63,rec=1) ((areas_mask(i,j,1),i=1,nx),j=1,ny)
!       read(63,rec=2) ((areas_mask(i,j,2),i=1,nx),j=1,ny)
      endif

! Define whether this simulation is going to restrict the
!   assimilation influence to some local area surrounding each
!   data point that is assimilated.  This was implemented for
!   ANWR simulations where we only had observations in a corner
!   of the simulation domain and I didn't want the corrections
!   to extend too far outside that local observation area.  So,
!   this is an example of what can be done, and it is not written
!   for general application.  If you want to do something similar,
!   the associated subroutine can be edited for your specific
!   simulation of interest.  For yes, local_assim_flag = 1, for
!   no, local_assim_flag = 0.
      local_assim_flag = 0

! Identify the file that contains the local data assimilation
!   mask.  This is only used if local_assim_flag = 1.
      fname_sweobs_barnes_mask = &
     &  '../swe_obs/2014/barnes/obs.gridded.gdat'

! END USER EDIT SECTION.
! END USER EDIT SECTION.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Open the swe input file.
      open (unit=61,file=fname_sweobs)

! Read in the number of years this assimilation will process.
      read(61,*) nyears

! LOOP THROUGH THE YEARS IN THE SIMULATION.
      do nyear=1,nyears

! Run the data assimilation routines.
        call data_assimilation(nx,ny,deltax,deltay,beta,&
     &    areas_flag,sprec_ratio,smelt_ratio,corr_factor,&
     &    areas_mask,xmn,ymn,xstn,ystn,nobs_dates,iobs_rec,&
     &    local_assim_flag,fname_swed,fname_sspr,fname_ssmt,&
     &    fname_sweobs,fname_sweobs_barnes_mask,iday_init,&
     &    imonth_init,iyear_init,nyear,nobs_total,print_inc,&
     &    nyears,dt,max_iter,xhour_init)

! Build an array indicating the appropriate correction factor to
!   use at any given time during the simulation.  What this does
!   is define an index array that contains the record number that
!   gets used at every model time step during the second model run
!   loop.  This record number corresponds to the record (krec) of
!   the corr_factor(i,j,krec) array that was generated and saved
!   in the subroutine above.
        call corr_factor_index(nobs_dates,icorr_factor_index,&
     &    iobs_rec,max_iter,sprec_ratio,smelt_ratio,print_inc,&
     &    iday_init,imonth_init,iyear_init,nyear,nobs_total_cfi,&
     &    nyears)

      enddo

      return
      end SUBROUTINE DATAASSIM_USER
!      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine data_assimilation(nx,ny,deltax,deltay,beta,&
     &  areas_flag,sprec_ratio,smelt_ratio,corr_factor,&
     &  areas_mask,xmn,ymn,xstn,ystn,nobs_dates,iobs_rec,&
     &  local_assim_flag,fname_swed,fname_sspr,fname_ssmt,&
     &  fname_sweobs,fname_sweobs_barnes_mask,iday_init,&
     &  imonth_init,iyear_init,nyear,nobs_total,print_inc,&
     &  nyears,dt,max_iter,xhour_init)

      use snowmodel_inc
      implicit none

      real deltax,deltay,undef,dn,beta,areas_flag,swe_count,cf_min,&
     &  print_inc,dt,xhour_init
      real sprec_ratio(max_obs_dates),smelt_ratio(max_obs_dates)
      real corr_factor(nx_max,ny_max,max_obs_dates)
      real swe_tmp(nx_max,ny_max),sum_sprec_tmp1(nx_max,ny_max),&
     &  sum_sprec_tmp2(nx_max,ny_max),grid(nx_max,ny_max),&
     &  areas_mask(nx_max,ny_max),sum_smelt_tmp1(nx_max,ny_max),&
     &  sum_smelt_tmp2(nx_max,ny_max)
      real corr_factor_tmp(nx_max*ny_max),swe_obs(nx_max*ny_max),&
     &  swe_model(nx_max*ny_max),sumsprec_model(nx_max*ny_max),&
     &  delta_old(nx_max*ny_max),obsid(nx_max*ny_max),&
     &  sumsmelt_model(nx_max*ny_max),obsid_old(nx_max*ny_max),&
     &  delta_old_tmp(nx_max*ny_max)
!     real corr_offset(nx_max*ny_max)

      double precision xmn,ymn
      double precision xstn(nx_max*ny_max),ystn(nx_max*ny_max)

      integer iobs_rec(max_obs_dates)
      integer ii(nx_max*ny_max),jj(nx_max*ny_max)
      integer iobs_num,irec1,irec2,nobs_dates,nx,ny,i,j,ifill,&
     &  iobsint,k,nstns,nstns_old,kk,local_assim_flag,iiyr,iimo,&
     &  iidy,iobs_rec_tmp,iday_init,imonth_init,iyear_init,nyear,&
     &  krec,nobs_total,nyears,max_iter

      character*89 fname_swed,fname_sspr,fname_ssmt
      character*80 fname_sweobs
      character*80 fname_sweobs_barnes_mask

! Perform some initialization steps.
      if (nyear.eq.1) then

! Define some of the constants and parameters used in the data
!   assimilation.  ifill should be = 1; in that case undef is not
!   really used (so it does not have to be the same as defined
!   in the .par file).
        undef = -9999.0
        ifill = 1
        iobsint = 0

! Open a file to write some basic correction factor information.
!   This just saves information that the user might want to look
!   at.
        open (unit=77,file='swe_assim/corr_factor.txt')

! Open an output file for the correction factor array.
        open(62,file='swe_assim/corr_factor.gdat',&
     &    form='unformatted',access='direct',recl=4*nx*ny)

! Initialize the number-of-observations counter.
        nobs_total = 0

      endif

! Read the number of observation dates for this year.
      read(61,*) nobs_dates

! If you have observations for this year, generate the correction
!   factors.  For the case of no observations for this year, set
!   the correction factor equal to 1.0.
      if (nobs_dates.gt.0) then

! Loop through the observation dates.
        do iobs_num=1,nobs_dates

! Increment the number-of-observations counter.
          nobs_total = nobs_total + 1

! Read the date corresponding to this observation.
          read(61,*) iiyr,iimo,iidy

! Convert this date to the corresponding record number in the
!   original SnowModel output files.  Note that this has assumed
!   that daily data files were written out since the start of
!   the simulation.
          call get_obs_record(iday_init,imonth_init,iyear_init,&
     &      iidy,iimo,iiyr,iobs_rec_tmp,print_inc,dt)
          iobs_rec(iobs_num) = iobs_rec_tmp

! For this observation date, read in the data describing the
!   location and swe values for each observation.  For areas
!   simulations, xstn, and ystn correspond to the center of the
!   area domain and they are not really used.
          read(61,*) nstns
          do k=1,nstns
            read(61,*) obsid(k),xstn(k),ystn(k),swe_obs(k)
          enddo

! Convert the x and y locations to (ii,jj) locations.
          do k=1,nstns
            ii(k) = 1 + nint((xstn(k) - xmn) / deltax)
            jj(k) = 1 + nint((ystn(k) - ymn) / deltay)
          enddo

! If you do a data assimilation run from start to finish, it is
!   not required to close and reopen these files.  But if you are
!   doing a history restart then these files are no longer open
!   so you must do this.  What I do below works for both cases.
          close (238)
          close (239)
          close (240)

! Open the required inputs from the initial assimilation loop.
!   Open swe depth (swe_depth).
!     /outputs/wo_assim/swed.gdat is unit 238 in outputs_user.f
!   Open sum snow precip (sum_sprec).
!     /outputs/wo_assim/sspr.gdat is unit 239 in outputs_user.f
!   Open sum snow melt (sum_smelt).
!     /outputs/wo_assim/ssmt.gdat is unit 240 in outputs_user.f
          open (238,file=fname_swed,&
     &      form='unformatted',access='direct',recl=4*1*nx*ny)
          open (239,file=fname_sspr,&
     &      form='unformatted',access='direct',recl=4*1*nx*ny)
          open (240,file=fname_ssmt,&
     &      form='unformatted',access='direct',recl=4*1*nx*ny)

! Read the model output for the first observation time.
          if (iobs_num.eq.1) then
            irec1 = iobs_rec(iobs_num)
            read(238,rec=irec1) ((swe_tmp(i,j),i=1,nx),j=1,ny)
            read(239,rec=irec1) ((sum_sprec_tmp1(i,j),i=1,nx),j=1,ny)
            read(240,rec=irec1) ((sum_smelt_tmp1(i,j),i=1,nx),j=1,ny)

! For points, just pull the data at the appropriate grid cell.
!   For areas, average the data over the masked out area for each
!   'station'.
            do k=1,nstns
              if (areas_flag.eq.0.0) then
                swe_model(k) = swe_tmp(ii(k),jj(k))
                sumsprec_model(k) = sum_sprec_tmp1(ii(k),jj(k))
                sumsmelt_model(k) = sum_smelt_tmp1(ii(k),jj(k))
              elseif (areas_flag.eq.1.0) then
                swe_model(k) = 0.0
                sumsprec_model(k) = 0.0
                sumsmelt_model(k) = 0.0
                swe_count = 0.0
                do j=1,ny
                  do i=1,nx
                    if (areas_mask(i,j).eq.obsid(k)) then
! The following is used if the mask changes with observation time.
!                   if (areas_mask(i,j,iobs_num).eq.obsid(k)) then
                      swe_count = swe_count + 1.0
                      swe_model(k) = swe_model(k) + swe_tmp(i,j)
                      sumsprec_model(k) = sumsprec_model(k) +&
     &                  sum_sprec_tmp1(i,j)
                      sumsmelt_model(k) = sumsmelt_model(k) +&
     &                  sum_smelt_tmp1(i,j)
                    endif
                  enddo
                enddo
                swe_model(k) = swe_model(k) / swe_count
                sumsprec_model(k) = sumsprec_model(k) / swe_count
                sumsmelt_model(k) = sumsmelt_model(k) / swe_count
              endif
            enddo
          endif

! Read the model output for any additional observation times (irec1
!   = current obs time, irec2 = previous obs time).
          if (iobs_num.gt.1) then
            irec1 = iobs_rec(iobs_num)
            irec2 = iobs_rec(iobs_num-1)
            read(238,rec=irec1) ((swe_tmp(i,j),i=1,nx),j=1,ny)
            read(239,rec=irec1) ((sum_sprec_tmp1(i,j),i=1,nx),j=1,ny)
            read(239,rec=irec2) ((sum_sprec_tmp2(i,j),i=1,nx),j=1,ny)
            read(240,rec=irec1) ((sum_smelt_tmp1(i,j),i=1,nx),j=1,ny)
            read(240,rec=irec2) ((sum_smelt_tmp2(i,j),i=1,nx),j=1,ny)

! For points, just pull the data at the appropriate grid cell.
!   For areas, average the data over the masked out area for each
!   'station'.
            do k=1,nstns
              if (areas_flag.eq.0.0) then
                swe_model(k) = swe_tmp(ii(k),jj(k))
                sumsprec_model(k) = sum_sprec_tmp1(ii(k),jj(k)) -&
     &            sum_sprec_tmp2(ii(k),jj(k))
                sumsmelt_model(k) = sum_smelt_tmp1(ii(k),jj(k)) -&
     &            sum_smelt_tmp2(ii(k),jj(k))
              elseif (areas_flag.eq.1.0) then
                swe_model(k) = 0.0
                sumsprec_model(k) = 0.0
                sumsmelt_model(k) = 0.0
                swe_count = 0.0
                do j=1,ny
                  do i=1,nx
                    if (areas_mask(i,j).eq.obsid(k)) then
! The following is used if the mask changes with observation time.
!                   if (areas_mask(i,j,iobs_num).eq.obsid(k)) then
                      swe_count = swe_count + 1.0
                      swe_model(k) = swe_model(k) + swe_tmp(i,j)
                      sumsprec_model(k) = sumsprec_model(k) + &
     &                  sum_sprec_tmp1(i,j) - sum_sprec_tmp2(i,j)
                      sumsmelt_model(k) = sumsmelt_model(k) + &
     &                  sum_smelt_tmp1(i,j) - sum_smelt_tmp2(i,j)
                    endif
                  enddo
                enddo
                swe_model(k) = swe_model(k) / swe_count
                sumsprec_model(k) = sumsprec_model(k) / swe_count
                sumsmelt_model(k) = sumsmelt_model(k) / swe_count
              endif
            enddo
          endif

! To avoid a divide by zero later on, make sure sumsprec_model and
!   sumsmelt_model are not both zero.
          do k=1,nstns
            sumsprec_model(k) = sumsprec_model(k) + 1.0e-6
          enddo

! Determine whether we will adjust the precipitation or melt.  To
!   do this, calculate the relative contributions of precipitation
!   and melt inputs for this correction period.  This can be
!   different for each observation interval.  Calculate the average
!   over all of the stations/areas in the domain.
          sprec_ratio(iobs_num) = 0.0
          smelt_ratio(iobs_num) = 0.0
          do k=1,nstns
            sprec_ratio(iobs_num) = sprec_ratio(iobs_num) + &
     &        sumsprec_model(k) / (sumsprec_model(k)+sumsmelt_model(k))
            smelt_ratio(iobs_num) = smelt_ratio(iobs_num) + &
     &        sumsmelt_model(k) / (sumsprec_model(k)+sumsmelt_model(k))
          enddo
          sprec_ratio(iobs_num) = sprec_ratio(iobs_num) / real(nstns)
          smelt_ratio(iobs_num) = smelt_ratio(iobs_num) / real(nstns)

! Initialize the delta swe variable.
          if (iobs_num.eq.1) then
            do k=1,nstns
              delta_old(k) = 0.0
            enddo
          else
            do k=1,nstns
              delta_old(k) = 0.0
            enddo
            do k=1,nstns
              do kk=1,nstns_old
                if(obsid(k).eq.obsid_old(kk)) &
     &            delta_old(k) = delta_old_tmp(kk)
              enddo
            enddo
!           write (77,*)
!           do k=1,nstns
!             write (77,*) 'k, delta_old(k)',k,100.*delta_old(k)
!           enddo
!           write (77,*)
          endif

! Calculate the correction factor to be used in the next model
!   iteration.  Let the correction factor equal 1.0 during
!   periods where we have no swe observations.  Also, note that the
!   reason for the delta_old variable is to account for the fact
!   that that delta will be fixed with the previous date correction
!   time period.  This is one of the things that allows the
!   correction to be done in two model iterations.
! If sumsprec_model or sumsmelt_model are too small to be used in
!   the assimilation (like less than 1 mm), set corr_factor_tmp = 1.0
!   so no adjustments are performed for this observation interval.
          cf_min = 0.1
          do k=1,nstns
            if (sprec_ratio(iobs_num).ge.smelt_ratio(iobs_num)) then
              if (sumsprec_model(k).lt.1.0e-3) then
                corr_factor_tmp(k) = 1.0
              else
                corr_factor_tmp(k) = 1.0 + &
     &            (swe_obs(k) - swe_model(k) - delta_old(k)) / &
     &            sumsprec_model(k)
                corr_factor_tmp(k) = max(cf_min,corr_factor_tmp(k))
              endif
            else
              if (sumsmelt_model(k).lt.1.0e-3) then
                corr_factor_tmp(k) = 1.0
              else
                corr_factor_tmp(k) = 1.0 + &
     &            (swe_model(k) - swe_obs(k) + delta_old(k)) / &
     &            sumsmelt_model(k)
                corr_factor_tmp(k) = max(cf_min,corr_factor_tmp(k))
              endif
            endif
! Save some information about the model calculations.
!           write (77,*) '---'
!           write (77,*) k,swe_obs(k)
!           write (77,*) k,swe_model(k)
!           write (77,*) k,delta_old(k)
!           write (77,*) k,swe_obs(k)-swe_model(k)-delta_old(k)
!           write (77,*) k,sumsprec_model(k)
!           write (77,*) k,sumsmelt_model(k)
!           write (77,*) k,corr_factor_tmp(k)
!           write (77,*) '---'

! Save some data from this observation time for use at the next
!   observation time.
            nstns_old = nstns
            obsid_old(k) = obsid(k)
            delta_old_tmp(k) = swe_obs(k) - swe_model(k)
          enddo

! Now that I have the correction factors calculated at each
!   observation point, interpolate those over the simulation domain.

! Use the barnes oi scheme to create the distribution. If there is
!   only a single station, distribute those data uniformly over
!   the domain.
          if (nstns.ge.2) then
            call get_dn(nx,ny,deltax,deltay,nstns,dn,iobsint)

! Modify the size of dn.
            dn = beta * dn

            call barnes_oi(nx,ny,deltax,deltay,xmn,ymn, &
     &        nstns,xstn,ystn,corr_factor_tmp,dn,grid,undef,ifill)
          elseif (nstns.eq.1) then
            call single_stn(nx,ny,nstns,corr_factor_tmp,grid)
          endif

! The following calculations are done if you want to implement the
!   special case where you limit the data assimilation corrections
!   to a certain area of your simulation domain.  Edits to this
!   subroutine will certainly be required to make it work for your
!   specific application.  This subroutine generates correction
!   factors with 1.0's in places outside the local obs influences.
          if (local_assim_flag.eq.1) then
            call mk_local_cfs(nx,ny,undef,xmn,ymn,deltax,deltay,&
     &        fname_sweobs,fname_sweobs_barnes_mask,nobs_dates,&
     &        corr_factor_tmp,beta,iobsint,ifill,grid)
          endif

! Define the correction surface record that corresponds to this
!   year and observation.
          krec = nobs_total + (nyear - 1)
          if (krec.gt.max_obs_dates) then
            print *, 'max_obs_dates must be increased in snowmodel.inc'
            print *, 'krec = ',krec,'  max_obs_dates = ',max_obs_dates
            stop
          endif

! Use the gridded output file to build the corr_factor array.
          do j=1,ny
            do i=1,nx
              corr_factor(i,j,krec) = grid(i,j)
              corr_factor(i,j,krec) = &
     &          max(cf_min,corr_factor(i,j,krec))
            enddo
          enddo

! Note that the interpolation scheme may have produced correction
!   factors that do not produce exact matches with the
!   observations (like happens with the case of having a single
!   data point).  If you are interested, calculate the difference
!   between the exact value and the actual calculated value, and
!   then write it out as done below.
!         do k=1,nstns
!           if (sprec_ratio(iobs_num).ge.smelt_ratio(iobs_num)) then
!             corr_offset(k) = sumsprec_model(k) *
!    &          (corr_factor(ii(k),jj(k),iobs_num) - corr_factor_tmp(k))
!           else
!             corr_offset(k) = sumsmelt_model(k) *
!    &          (corr_factor(ii(k),jj(k),iobs_num) - corr_factor_tmp(k))
!           endif
!         enddo

! Write some information to the text file.
          write (77,*) &
     &  '**************************************************************'

!         write (77,*) ' sprec_ratio =',sprec_ratio(iobs_num),
!    &      '  smelt_ratio =',smelt_ratio(iobs_num)
!         write (77,*)

          write (77,*) iiyr,iimo,iidy
          write (77,*) '          sprec_ratio =',sprec_ratio(iobs_num)
          write (77,*) '          smelt_ratio =',smelt_ratio(iobs_num)
          write (77,*)

          do k=1,nstns
!           write (77,*) k,' swe diff =',
!    &        100.0*abs(swe_obs(k)-swe_model(k)),' SWE OBS =',
!    &        100.0*swe_obs(k)
!           write (77,*) 'sumsprec =',sumsprec_model(k)*100.,
!    &        '  SWE MODEL =',swe_model(k)*100.
!           write (77,*) 'iobs_num =',iobs_num,
!    &        '  CORR_FACTOR =',corr_factor_tmp(k)

            write (77,*) '         SWE OBS (cm) =',100.0*swe_obs(k)
            write (77,*) '       SWE MODEL (cm) =',100.0*swe_model(k)
            write (77,*) '          CORR_FACTOR =',corr_factor_tmp(k)

!c          write (77,*) 'corr_offset =',100.*corr_offset(k),
!c   &        '  ij',ii(k),jj(k)
!c          write (77,*) '     delta_old =',100.*delta_old(k),
!c   &        '      corr fact used =',corr_factor(ii(k),jj(k),iobs_num)
!           write (77,*)
!           write (77,*) k,' sumsprec_model(k) =',sumsprec_model(k)
!           write (77,*) k,' sumsmelt_model(k) =',sumsmelt_model(k)
!           write (77,*)
          enddo

          write (77,*) &
     &  '**************************************************************'

! Write the output data to a grads file.
          write(62,rec=krec) ((corr_factor(i,j,krec),i=1,nx),j=1,ny)

        enddo

! Fill corr_factor with 1.0 for the period following the last obs
!   date in the current year.  This is also required for the history
!   restart to work correctly.  Without the history restart this was
!   already done as part of the model initialization.
        if (krec+1.gt.max_obs_dates) then
          print *, 'max_obs_dates must be increased in snowmodel.inc'
          print *, 'krec+1 = ',krec+1,'  max_obs_dates = ',max_obs_dates
          stop
        endif
        do j=1,ny
          do i=1,nx
            corr_factor(i,j,krec+1) = 1.0
          enddo
        enddo

        write(62,rec=krec+1) ((corr_factor(i,j,krec+1),i=1,nx),j=1,ny)

! The met, topo, and veg files must be closed for the next model
!   iteration.
        close (20)
        close (37)
        close (38)

        close (238)
        close (239)
        close (240)

      else

! For the case of no observations for this year, set the correction
!   factor equal to 1.0.
        krec = nobs_total + nyear
        if (krec.gt.max_obs_dates) then
          print *, 'max_obs_dates must be increased in snowmodel.inc'
          print *, 'krec = ',krec,'  max_obs_dates = ',max_obs_dates
          stop
        endif

        do j=1,ny
          do i=1,nx
            corr_factor(i,j,krec) = 1.0
          enddo
        enddo

        write(62,rec=krec) ((corr_factor(i,j,krec),i=1,nx),j=1,ny)

      endif

! Create the GrADS .ctl (control) file to go with the GrADS
!   .gdat corr_factor file that was generated by this model run.
      call mk_cf_prec_ctl(nx,ny,deltax,deltay,xmn,ymn,dt, &
     &  iyear_init,imonth_init,iday_init,xhour_init,max_iter, &
     &  nobs_total,nyears)

      return
      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine corr_factor_index(nobs_dates,icorr_factor_index, &
     &  iobs_rec,max_iter,sprec_ratio,smelt_ratio,print_inc, &
     &  iday_init,imonth_init,iyear_init,nyear,nobs_total_cfi, &
     &  nyears)

      use snowmodel_inc
      implicit none

      integer icorr_factor_index(max_time_steps)

      integer kk,istart,iend,nobs_dates,iter,max_iter,krec,nyear,&
     &  nobs_total_cfi,ioptn,julian_start,iiyr,julian_end,iday_init,&
     &  imonth_init,iyear_init,nyears
      integer iobs_rec(max_obs_dates)
      real sprec_ratio(max_obs_dates),smelt_ratio(max_obs_dates) 
      real print_inc

! Initialize the number-of-observations counter.
      if (nyear.eq.1) nobs_total_cfi = 0

! Build an array indicating the appropriate correction factor to
!   use at each time step during the simulation.
      if (nobs_dates.gt.0) then

! Loop through the observation dates.
        do kk=1,nobs_dates+1

! Increment the number-of-observations counter.
          if (kk.le.nobs_dates) nobs_total_cfi = nobs_total_cfi + 1

! FIRST, FROM THE SIMULATION START UNTIL THE FIRST OBSERVATION, FOR
!   EACH YEAR.
          if (kk.eq.1) then

! Here istart equals the first model time step of each year.
            ioptn = 3
            call calndr (ioptn,iday_init,imonth_init,iyear_init,&
     &        julian_start)

! Find the Julian day for this data record.
            iiyr = iyear_init + (nyear - 1)
            call calndr (ioptn,iday_init,imonth_init,iiyr,julian_end)

! Calculate istart.
            istart = 1 + (julian_end - julian_start) * nint(print_inc)

! Here iend equals the model time step corresponding to the end of
!   the first observation day.  Take advantage of the data-file
!   record that was calculated before.  Convert from daily data
!   output records to model time steps using print_inc.
            iend = iobs_rec(kk) * nint(print_inc)

! Define the corr_factor data array record.
            krec = nobs_total_cfi + (nyear - 1)

! Fill the index for each model time step.
            do iter=istart,iend
              if (sprec_ratio(kk).ge.smelt_ratio(kk)) then
                icorr_factor_index(iter) = krec
              else
                icorr_factor_index(iter) = -krec
              endif
            enddo

! SECOND, BETWEEN THE LAST OBSERVATION AND THE END OF THE SIMULATION.
          elseif (kk.eq.nobs_dates+1) then
            istart = iobs_rec(kk-1) * nint(print_inc) + 1

! Here iend equals the last time step in the year of interest.
! Find the Julian day for this data record.
            iiyr = iyear_init + nyear
            call calndr (ioptn,iday_init,imonth_init,iiyr,julian_end)

! Calculate iend for this year.
            iend = (julian_end - julian_start) * nint(print_inc)
            iend = min(iend,max_iter)

! Define the corr_factor data array record.
            krec = nobs_total_cfi + nyear

! Fill the index for each model time step.
            do iter=istart,iend
              icorr_factor_index(iter) = krec
            enddo

! THIRD, ANY PERIODS BETWEEN OBSERVATIONS.
          else
            istart = iobs_rec(kk-1) * nint(print_inc) + 1
            iend = iobs_rec(kk) * nint(print_inc)

! Define the corr_factor data array record.
            krec = nobs_total_cfi + (nyear - 1)

! Fill the index for each model time step.
            do iter=istart,iend
              if (sprec_ratio(kk).ge.smelt_ratio(kk)) then
                icorr_factor_index(iter) = krec
              else
                icorr_factor_index(iter) = -krec
              endif
            enddo
          endif
        enddo

      else

! Create an array indes for the case of no observations for this
!   year.  Here istart equals the first model time step of each year.
        ioptn = 3
        call calndr (ioptn,iday_init,imonth_init,iyear_init,&
     &    julian_start)

! Find the Julian day for this data record.
        iiyr = iyear_init + (nyear - 1)
        call calndr (ioptn,iday_init,imonth_init,iiyr,julian_end)

! Calculate istart.
        istart = 1 + (julian_end - julian_start) * nint(print_inc)

! Calculate iend for this year.
        iiyr = iyear_init + nyear
        call calndr (ioptn,iday_init,imonth_init,iiyr,julian_end)
        iend = (julian_end - julian_start) * nint(print_inc)
        iend = min(iend,max_iter)

! Define the corr_factor data array record.
        krec = nobs_total_cfi + nyear

! Fill the index for each model time step.
        do iter=istart,iend
          icorr_factor_index(iter) = krec
        enddo

      endif

! SAVE A VERSION OF THE INDEX THAT THE USER CAN EASILY LOOK AT.
      if (nyear.eq.nyears) then
        print *
        do iter=1,max_iter
          write (77,*) iter,icorr_factor_index(iter)
        enddo
        print *
      endif

      return
      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine mk_local_cfs(nx,ny,undef,xmn,ymn,deltax,deltay,&
     &  fname_sweobs,fname_sweobs_barnes_mask,nobs_dates,&
     &  corr_factor_tmp,beta,iobsint,ifill,grid)

      use snowmodel_inc
      implicit none

      integer nstns,k,nx,ny,i,j,ntstations,icount,nobs_dates,ifill, &
     &  iobsint,nstns2

      real dummy,stnid,undef,deltax,deltay,beta,dn
      real xmask(nx_max,ny_max),grid(nx_max,ny_max)
      real corr_factor_tmp(nx_max*ny_max),cf_tmp2(nx_max*ny_max),&
     &  obsid(nx_max*ny_max)

      real var(nstns_max)

      double precision xmn,ymn
      double precision yg(nx_max,ny_max),xg(nx_max,ny_max)
      double precision xstn(nstns_max),ystn(nstns_max)

      character*80 fname_sweobs
      character*80 fname_sweobs_barnes_mask

      print *
      print *,'You are doing local assimilation, this requires'
      print *,'an observation mask to have been generated before'
      print *,'the model run and saved in the file called:'
      print *,'fname_sweobs_barnes_mask'
      print *
      print *,'This was also not tested after some big changes to'
      print *,'the data assimilation code.  So I strongly suggest'
      print *,'you make sure this is doing what you want when you'
      print *,'first start using it.'
      print *
      stop
      if (nobs_dates.gt.1) then
        print *,'THIS HAS NOT BEEN MADE TO WORK WITH MORE THAN'
        print *,'ONE OBS TIME.'
        stop
      endif
      print *

! Save the correction factors for the local-influence assimilation
!   scheme.
      open (unit=78,file='swe_assim/cf.txt')
      do k=1,nstns
        write (78,*) corr_factor_tmp(k)
      enddo
      close (78)

! These are the obs and correction factors calculated from the
!   first loop in SnowModel.
      open (721,file=fname_sweobs)
      open (731,file='swe_assim/cf.txt')

      read (721,*) nstns
      do k=1,nstns
        read (721,*) stnid,xstn(k),ystn(k),dummy
        read (731,*) var(k)
!       if (var(k).gt.1.5) then
!         print *, 'cf > 1.5 found; setting to 1.5',k,var(k)
!         var(k) = 1.5
!       endif
!       if (var(k).lt.0.5) then
!         print *, 'cf < 0.5 found; setting to 0.5',k,var(k)
!         var(k) = 0.5
!       endif
!       var(k) = min(1.5,var(k))
!       var(k) = max(0.5,var(k))
      enddo
      close (721)
      close (731)

! Create a collection of 'stations' with correction factors of
!   1.0 in areas outside of our traverse regions.
      open(741,file=fname_sweobs_barnes_mask,&
     &  form='unformatted',access='direct',recl=4*nx*ny)
      read (741,rec=1) ((xmask(i,j),i=1,nx),j=1,ny)
      close (741)

! Create an array of e, n coordinates for this domain.
      do j=1,ny
        do i=1,nx
          xg(i,j) = xmn + deltax * (real(i) - 1.0)
          yg(i,j) = ymn + deltay * (real(j) - 1.0)
        enddo
      enddo

      do j=1,ny
        do i=1,nx
          if (xmask(i,j).ne.undef) then
            xg(i,j) = undef
            yg(i,j) = undef
          endif
        enddo
      enddo

! Count how many cf=1.0 'stations' you are going to end up with.
      icount = 0
      do j=1,ny,100
        do i=1,nx,100
          if (xg(i,j).ne.undef) then
            icount = icount + 1
          endif
        enddo
      enddo

! Write out the original stations.
      open (761,file='swe_assim/cf_with_mask.txt')
      ntstations = nstns + icount
      write (761,88) ntstations
      do k=1,nstns
        write (761,89) k,xstn(k),ystn(k),var(k)
      enddo

! Write out the cf=1.0 stations.
      icount = 0
      do j=1,ny,100
        do i=1,nx,100
          if (xg(i,j).ne.undef) then
            icount = icount + 1
            write (761,89) icount+1000,xg(i,j),yg(i,j),1.0
          endif
        enddo
      enddo
      close (761)

! Read in the new local cf data.
      open (79,file='swe_assim/cf_with_mask.txt')
      read (79,*) nstns2
      do k=1,nstns2
        read (79,*) obsid(k),xstn(k),ystn(k),cf_tmp2(k)
      enddo
      close (79)

      call get_dn(nx,ny,deltax,deltay,nstns2,dn,iobsint)

      dn = beta * dn

      call barnes_oi(nx,ny,deltax,deltay,xmn,ymn,&
     &  nstns2,xstn,ystn,cf_tmp2,dn,grid,undef,ifill)

! Write the output data to a grads file.
      open(511,file='swe_assim/corr_factor_w-mask.gdat',&
     &  form='unformatted',access='direct',recl=4*nx*ny)
      write(511,rec=1) ((grid(i,j),i=1,nx),j=1,ny)
      close (511)

  88  format (i10)
  89  format (i10,2f15.1,f10.4)

      return
      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine get_obs_record(iday_init,imonth_init,iyear_init,&
     &  iidy,iimo,iiyr,iobs_rec_tmp,print_inc,dt)

      implicit none

      integer ioptn,iday_init,imonth_init,iyear_init,julian_start,&
     &  iidy,iimo,iiyr,julian_end,iobs_rec_tmp,n_writes_per_day

      real print_inc,dt

! Find the Julian day at the start of the model run.
      ioptn = 3
      call calndr (ioptn,iday_init,imonth_init,iyear_init,julian_start)

! Find the Julian day for this data record.
      call calndr (ioptn,iidy,iimo,iiyr,julian_end)

! Calculate the day of simulation for this data record.  This is the
!   same as the output file data record.
      iobs_rec_tmp = julian_end - julian_start + 1

! Correct this for the case where sub-daily data writes were made.
      n_writes_per_day = nint((86400.0 / dt) / print_inc)
      iobs_rec_tmp = n_writes_per_day * iobs_rec_tmp

      return
      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine mk_cf_prec_ctl(nx,ny,deltax,deltay,xmn,ymn,dt,&
     &  iyear_init,imonth_init,iday_init,xhour_init,max_iter,&
     &  nobs_total,nyears)

      implicit none

      integer nx,ny,iyear_init,imonth_init,iday_init,max_iter,&
     &  igrads_dt,nobs_total,nyears
      double precision xmn,ymn,xmn_km,ymn_km
      real deltax,deltay,deltax_km,deltay_km,xhour_init,dt,undef

      integer len_desc,trailing_blanks

      character*95 output_fname
      character*80 filename,description
      character*3 cmo(12)
      character*2 cdt

      data cmo /'jan','feb','mar','apr','may','jun',&
     &          'jul','aug','sep','oct','nov','dec'/

      data description /&
     &  'cf    0  0 precip correction factor (above and below 1.0)'/

      undef = -9999.0

      deltax_km = deltax / 1000.0
      deltay_km = deltay / 1000.0
      xmn_km = xmn / 1000.0
      ymn_km = ymn / 1000.0

! Convert the write interval from seconds to hours or a day.
      igrads_dt = 1
      cdt = 'dy'

      filename = 'swe_assim/corr_factor.ctl'
      output_fname = 'corr_factor.gdat'

      len_desc = 80 - trailing_blanks(description)

      open (71,file=filename)

      write (71,51) output_fname

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
      write (71,57) nobs_total+nyears,nint(xhour_init),&
     &  iday_init,cmo(imonth_init),iyear_init,igrads_dt,cdt
      write (71,58)
      write (71,59) description(1:len_desc)
      write (71,60)

      close (71)

! This "a" by itself clips the trailing blanks in the a80 string.
   51 format ('DSET ^',a)
   52 format ('TITLE SnowModel data assimilation precip corr factor')
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! file = calendar.f  version 1.0
!
! Note from Glen: This is the version that should be used as
!   part of my computer programs!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This program performs various date conversion calculations
! using subroutine calndr().
! On a 32-bit computer, it can handle any date between roughly
! 5 million BC and 5 million AD.  This limitation is due to
! the range of integers that can be expressed with 32 bits.
! The algorithm has no limitation.
!
! Using function idaywk(), the day of the week is computed
! along with the answer to the user's calendar calculation.
!
! External routines called:
! calndr  calendar conversions
! idaywk  day of the week determination
!
! Portability
! This routine is coded to Fortran 77 standards except that
! lower case is used.
!
! Copyright (C) 1999 Jon Ahlquist.
! Issued under the second GNU General Public License.
! See www.gnu.org for details.
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
! If you find any errors, please notify:
! Jon Ahlquist <ahlquist@met.fsu.edu>   
! Dept of Meteorology
! Florida State University
! Tallahassee, FL 32306-4520
! 15 March 1999.
!
!----------
!
! Declare variables.
!     implicit     none
!     integer      day, daynum, ioptn, julian, month, year
!     character*9  daynam(0:6)
!
! Declare the integer function used to compute the day of the week.
!     integer  idaywk
!
! Define the day names.
!     data  daynam /'Sunday',   'Monday', 'Tuesday', 'Wednesday',
!    &              'Thursday', 'Friday', 'Saturday'/
!
! Variables and their meanings
! day     day of the month.
! daynam  array of day names.  (daynam(0)='Sunday', daynam(1)='Monday',
!            ..., daynam(6)='Saturday')
! daynum  day number during the year.  (1 for 1 January, 2 for
!            2 January, 32 for 1 February, etc.)
! idaywk  integer function that returns an integer counter indicating
!            the day of the week, where 0 refers to Sunday, 1 to Monday,
!            up to 6 for Saturday.
! ioptn   option indicator where 0 < abs(ioptn) < 6.
!            See below and especially subroutine calndr for details.
! julian  Julian Day number.
! month   month counter (1=January, 2=February, ..., 12=December)
! year    year expressed with ALL digits.  DO NOT abbreviate years
!            by using only the last two digits.
!
!----------
!
      subroutine calndr (ioptn, iday, month, iyear, idayct)
!
!----------
!
! CALNDR = CALeNDaR conversions, version 1.0
!
! Input variable specifying the desired calendar conversion option.
      integer ioptn
!
! Input/Output variables (sometimes input, sometimes output,
! depending on the value of the desired option, ioptn.)
      integer  iday, month, iyear, idayct
!
!----------
!
! Subroutine calndr() performs calendar calculations using either
! the standard Gregorian calendar or the old Julian calendar.
! This subroutine extends the definitions of these calendar systems
! to any arbitrary year.  The algorithms in this subroutine
! will work with any date in the past or future,
! but overflows will occur if the numbers are sufficiently large.
! For a computer using a 32-bit integer, this routine can handle
! any date between roughly 5.8 million BC and 5.8 million AD
! without experiencing overflow during calculations.
!
! No external functions or subroutines are called.
!
!----------
!
! INPUT/OUTPUT ARGUMENTS FOR SUBROUTINE CALNDR()
!
! "ioptn" is the desired calendar conversion option explained below.
! Positive option values use the standard modern Gregorian calendar.
! Negative option values use the old Julian calendar which was the
! standard in Europe from its institution by Julius Caesar in 45 BC
! until at least 4 October 1582.  The Gregorian and Julian calendars
! are explained further below.
!
! (iday,month,iyear) is a calendar date where "iday" is the day of
! the month, "month" is 1 for January, 2 for February, etc.,
! and "iyear" is the year.  If the year is 1968 AD, enter iyear=1968,
! since iyear=68 would refer to 68 AD.
! For BC years, iyear should be negative, so 45 BC would be iyear=-45.
! By convention, there is no year 0 under the BC/AD year numbering
! scheme.  That is, years proceed as 2 BC, 1 BC, 1 AD, 2 AD, etc.,
! without including 0.  Subroutine calndr() will print an error message
! and stop if you specify iyear=0.
!
! "idayct" is a day count.  It is either the day number during the
! specified year or the Julian Day number, depending on the value
! of ioptn.  By day number during the specified year, we mean
! idayct=1 on 1 January, idayct=32 on 1 February, etc., to idayct=365
! or 366 on 31 December, depending on whether the specified year
! is a leap year.
!
! The values of input variables are not changed by this subroutine.
!
!
! ALLOWABLE VALUES FOR "IOPTN" and the conversions they invoke.
! Positive option values ( 1 to  5) use the standard Gregorian calendar.
! Negative option values (-1 to -5) use the old      Julian    calendar.
!
! Absolute
!  value
! of ioptn   Input variable(s)     Output variable(s)
!
!    1       iday,month,iyear      idayct
! Given a calendar date (iday,month,iyear), compute the day number
! (idayct) during the year, where 1 January is day number 1 and
! 31 December is day number 365 or 366, depending on whether it is
! a leap year.
!
!    2       idayct,iyear          iday,month
! Given the day number of the year (idayct) and the year (iyear),
! compute the day of the month (iday) and the month (month).
!
!    3       iday,month,iyear      idayct
! Given a calendar date (iday,month,iyear), compute the Julian Day
! number (idayct) that starts at noon of the calendar date specified.
!
!    4       idayct                iday,month,iyear
! Given the Julian Day number (idayct) that starts at noon,
! compute the corresponding calendar date (iday,month,iyear).
!
!    5       idayct                iday,month,iyear
! Given the Julian Day number (idayct) that starts at noon,
! compute the corresponding day number for the year (iday)
! and year (iyear).  On return from calndr(), "month" will always
! be set equal to 1 when ioptn=5.
!
! No inverse function is needed for ioptn=5 because it is
! available through option 3.  One simply calls calndr() with:
! ioptn = 3,
! iday  = day number of the year instead of day of the month,
! month = 1, and
! iyear = whatever the desired year is.
!
!----------
!
! EXAMPLES
! The first 6 examples are for the standard Gregorian calendar.
! All the examples deal with 15 October 1582, which was the first day
! of the Gregorian calendar.  15 October is the 288-th day of the year.
! Julian Day number 2299161 began at noon on 15 October 1582.
!
! Find the day number during the year on 15 October 1582
!     ioptn = 1
!     call calndr (ioptn, 15, 10, 1582,  idayct)
! calndr() should return idayct=288
!
! Find the day of the month and month for day 288 in year 1582.
!     ioptn = 2
!     call calndr (ioptn, iday, month, 1582, 288)
! calndr() should return iday=15 and month=10.
!
! Find the Julian Day number for 15 October 1582.
!     ioptn = 3
!     call calndr (ioptn, 15, 10, 1582, julian)
! calndr() should return julian=2299161
!
! Find the Julian Day number for day 288 during 1582 AD.
! When the input is day number of the year, one should specify month=1
!     ioptn = 3
!     call calndr (ioptn, 288, 1, 1582, julian)
! calndr() should return dayct=2299161
!
! Find the date for Julian Day number 2299161.
!     ioptn = 4
!     call calndr (ioptn, iday, month, iyear, 2299161)
! calndr() should return iday=15, month=10, and iyear=1582
! 
! Find the day number during the year (iday) and year
! for Julian Day number 2299161.
!     ioptn = 5
!     call calndr (ioptn, iday, month, iyear, 2299161)
! calndr() should return iday=288, month=1, iyear=1582
!
! Given 15 October 1582 under the Gregorian calendar,
! find the date (idayJ,imonthJ,iyearJ) under the Julian calendar.
! To do this, we call calndr() twice, using the Julian Day number
! as the intermediate value.
!     call calndr ( 3, 15,        10, 1582,    julian)
!     call calndr (-4, idayJ, monthJ, iyearJ,  julian)
! The first call to calndr() should return julian=2299161, and
! the second should return idayJ=5, monthJ=10, iyearJ=1582
!
!----------
!
! BASIC CALENDAR INFORMATION
!
! The Julian calendar was instituted by Julius Caesar in 45 BC.
! Every fourth year is a leap year in which February has 29 days.
! That is, the Julian calendar assumes that the year is exactly
! 365.25 days long.  Actually, the year is not quite this long.
! The modern Gregorian calendar remedies this by omitting leap years
! in years divisible by 100 except when the year is divisible by 400.
! Thus, 1700, 1800, and 1900 are leap years under the Julian calendar
! but not under the Gregorian calendar.  The years 1600 and 2000 are
! leap years under both the Julian and the Gregorian calendars.
! Other years divisible by 4 are leap years under both calendars,
! such as 1992, 1996, 2004, 2008, 2012, etc.  For BC years, we recall
! that year 0 was omitted, so 1 BC, 5 BC, 9 BC, 13 BC, etc., and 401 BC,
! 801 BC, 1201 BC, etc., are leap years under both calendars, while
! 101 BC, 201 BC, 301 BC, 501 BC, 601 BC, 701 BC, 901 BC, 1001 BC,
! 1101 BC, etc., are leap years under the Julian calendar but not
! the Gregorian calendar.
!
! The Gregorian calendar is named after Pope Gregory XIII.  He declared
! that the last day of the old Julian calendar would be Thursday,
! 4 October 1582 and that the following day, Friday, would be reckoned
! under the new calendar as 15 October 1582.  The jump of 10 days was
! included to make 21 March closer to the spring equinox.
!
! Only a few Catholic countries (Italy, Poland, Portugal, and Spain)
! switched to the Gregorian calendar on the day after 4 October 1582.
! It took other countries months to centuries to change to the
! Gregorian calendar.  For example, England's first day under the
! Gregorian calendar was 14 September 1752.  The same date applied to
! the entire British empire, including America.  Japan, Russia, and many
! eastern European countries did not change to the Gregorian calendar
! until the 20th century.  The last country to change was Turkey,
! which began using the Gregorian calendar on 1 January 1927.
!
! Therefore, between the years 1582 and 1926 AD, you must know
! the country in which an event was dated to interpret the date
! correctly.  In Sweden, there was even a year (1712) when February
! had 30 days.  Consult a book on calendars for more details
! about when various countries changed their calendars.
!
! DAY NUMBER DURING THE YEAR
! The day number during the year is simply a counter equal to 1 on
! 1 January, 32 on 1 February, etc., thorugh 365 or 366 on 31 December,
! depending on whether the year is a leap year.  Sometimes this is
! called the Julian Day, but that term is better reserved for the
! day counter explained below.
!
! JULIAN DAY NUMBER
! The Julian Day numbering system was designed by Joseph Scaliger
! in 1582 to remove ambiguity caused by varying calendar systems.
! The name "Julian Day" was chosen to honor Scaliger's father,
! Julius Caesar Scaliger (1484-1558), an Italian scholar and physician
! who lived in France.  Because Julian Day numbering was especially
! designed for astronomers, Julian Days begin at noon so that the day
! counter does not change in the middle of an astronmer's observing
! period.  Julian Day 0 began at noon on 1 January 4713 BC under the
! Julian calendar.  A modern reference point is that 23 May 1968
! (Gregorian calendar) was Julian Day 2,440,000.
!
! JULIAN DAY NUMBER EXAMPLES
!
! The table below shows a few Julian Day numbers and their corresponding
! dates, depending on which calendar is used.  A negative 'iyear' refers
! to BC (Before Christ).
!
!                     Julian Day under calendar:
! iday  month   iyear     Gregorian   Julian
!  24     11   -4714            0        -38
!   1      1   -4713           38          0
!   1      1       1      1721426    1721424
!   4     10    1582      2299150    2299160
!  15     10    1582      2299161    2299171
!   1      3    1600      2305508    2305518
!  23      5    1968      2440000    2440013
!   5      7    1998      2451000    2451013
!   1      3    2000      2451605    2451618
!   1      1    2001      2451911    2451924
!
! From this table, we can see that the 10 day difference between the
! two calendars in 1582 grew to 13 days by 1 March 1900, since 1900 was
! a leap year under the Julian calendar but not under the Gregorian
! calendar.  The gap will widen to 14 days after 1 March 2100 for the
! same reason.
! 
!----------
!
! PORTABILITY
!
! This subroutine is written in standard FORTRAN 77.
! It calls no external functions or subroutines and should run
! without problem on any computer having a 32-bit word or longer.
! 
!----------
!
! ALGORITHM
!
! The goal in coding calndr() was clear, clean code, not efficiency.
! Calendar calculations usually take a trivial fraction of the time
! in any program in which dates conversions are involved.
! Data analysis usually takes the most time.
!
! Standard algorithms are followed in this subroutine.  Internal to
! this subroutine, we use a year counter "jyear" such that
!  jyear=iyear   when iyear is positive
!       =iyear+1 when iyear is negative.
! Thus, jyear does not experience a 1 year jump like iyear does
! when going from BC to AD.  Specifically, jyear=0 when iyear=-1,
! i.e., when the year is 1 BC.
!
! For simplicity in dealing with February, inside this subroutine,
! we let the year begin on 1 March so that the adjustable month,
! February is the last month of the year.
! It is clear that the calendar used to work this way because the
! months September, October, November, and December refer to
! 7, 8, 9, and 10.  For consistency, jyear is incremented on 1 March
! rather than on 1 January.  Of course, everything is adjusted back to
! standard practice of years beginning on 1 January before answers
! are returned to the routine that calls calndr().
!
! Lastly, we use a trick to calculate the number of days from 1 March
! until the end of the month that precedes the specified month.
! That number of days is int(30.6001*(month+1))-122,
! where 30.6001 is used to avoid the possibility of round-off and
! truncation error.  For example, if 30.6 were used instead,
! 30.6*5 should be 153, but round-off error could make it 152.99999,
! which would then truncated to 152, causing an error of 1 day.
!
! Algorithm reference:
! Dershowitz, Nachum and Edward M. Reingold, 1990: Calendrical
! Calculations.  Software-Practice and Experience, vol. 20, number 9
! (September 1990), pp. 899-928.
!
! Copyright (C) 1999 Jon Ahlquist.
! Issued under the second GNU General Public License.
! See www.gnu.org for details.
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
! If you find any errors, please notify:
! Jon Ahlquist <ahlquist@met.fsu.edu>   
! Dept of Meteorology
! Florida State University
! Tallahassee, FL 32306-4520
! 15 March 1999.
!
!-----
! Declare internal variables.
      integer  jdref,  jmonth, jyear, leap,&
     &         n1yr, n4yr, n100yr, n400yr,&
     &         ndays, ndy400, ndy100, nyrs,&
     &         yr400, yrref
!
! Explanation of all internal variables.
! jdref   Julian Day on which 1 March begins in the reference year.
! jmonth  Month counter which equals month+1 if month .gt. 2
!          or month+13 if month .le. 2.
! jyear   Year index,  jyear=iyear if iyear .gt. 0, jyear=iyear+1
!            if iyear .lt. 0.  Thus, jyear does not skip year 0
!            like iyear does between BC and AD years.
! leap    =1 if the year is a leap year, =0 if not.
! n1yr    Number of complete individual years between iyear and
!            the reference year after all 4, 100,
!            and 400 year periods have been removed.
! n4yr    Number of complete 4 year cycles between iyear and
!            the reference year after all 100 and 400 year periods
!            have been removed.
! n100yr  Number of complete 100 year periods between iyear and
!            the reference year after all 400 year periods
!            have been removed.
! n400yr  Number of complete 400 year periods between iyear and
!            the reference year.
! ndays   Number of days since 1 March during iyear.  (In intermediate
!            steps, it holds other day counts as well.)
! ndy400  Number of days in 400 years.  Under the Gregorian calendar,
!            this is 400*365 + 100 - 3 = 146097.  Under the Julian
!            calendar, this is 400*365 + 100 = 146100.
! ndy100  Number of days in 100 years,  Under the Gregorian calendar,
!            this is 100*365 + 24 = 36524.   Under the Julian calendar,
!            this is 100*365 + 25 = 36525.
! nyrs    Number of years from the beginning of yr400
!              to the beginning of jyear.  (Used for option +/-3).
! yr400   The largest multiple of 400 years that is .le. jyear.
!
!
!----------------------------------------------------------------
! Do preparation work.
!
! Look for out-of-range option values.
      if ((ioptn .eq. 0) .or. (abs(ioptn) .ge. 6)) then
         write(*,*)'For calndr(), you specified ioptn = ', ioptn
         write(*,*) &
     &   'Allowable values are 1 to 5 for the Gregorian calendar'
         write(*,*) &
     &   'and -1 to -5 for the Julian calendar.'
         stop
      endif
!
! Options 1-3 have "iyear" as an input value.
! Internally, we use variable "jyear" that does not have a jump
! from -1 (for 1 BC) to +1 (for 1 AD).
      if (abs(ioptn) .le. 3) then
         if (iyear .gt. 0) then
            jyear = iyear
         elseif (iyear .eq. 0) then
            write(*,*) &
     &      'For calndr(), you specified the nonexistent year 0'
            stop
         else
            jyear = iyear + 1
         endif
!
!        Set "leap" equal to 0 if "jyear" is not a leap year
!        and equal to 1 if it is a leap year.
         leap = 0
         if ((jyear/4)*4 .eq. jyear) then
            leap = 1
         endif
         if ((ioptn .gt. 0)               .and. &
     &       ((jyear/100)*100 .eq. jyear) .and. &
     &       ((jyear/400)*400 .ne. jyear)      ) then
               leap = 0
         endif
      endif
!
! Options 3-5 involve Julian Day numbers, which need a reference year
! and the Julian Days that began at noon on 1 March of the reference
! year under the Gregorian and Julian calendars.  Any year for which
! "jyear" is divisible by 400 can be used as a reference year.
! We chose 1600 AD as the reference year because it is the closest
! multiple of 400 to the institution of the Gregorian calendar, making
! it relatively easy to compute the Julian Day for 1 March 1600
! given that, on 15 October 1582 under the Gregorian calendar,
! the Julian Day was 2299161.  Similarly, we need to do the same
! calculation for the Julian calendar.  We can compute this Julian
! Day knwoing that on 4 October 1582 under the Julian calendar,
! the Julian Day number was 2299160.  The details of these calculations
! is next. 
!    From 15 October until 1 March, the number of days is the remainder
! of October plus the days in November, December, January, and February:
! 17+30+31+31+28 = 137, so 1 March 1583 under the Gregorian calendar
! was Julian Day 2,299,298.  Because of the 10 day jump ahead at the
! switch from the Julian calendar to the Gregorian calendar, 1 March
! 1583 under the Julian calendar was Julian Day 2,299,308.  Making use
! of the rules for the two calendar systems, 1 March 1600 was Julian
! Day 2,299,298 + (1600-1583)*365 + 5 (due to leap years) =
! 2,305,508 under the Gregorian calendar and day 2,305,518 under the
! Julian calendar.
!    We also set the number of days in 400 years and 100 years.
! For reference, 400 years is 146097 days under the Gregorian calendar
! and 146100 days under the Julian calendar.  100 years is 36524 days
! under the Gregorian calendar and 36525 days under the Julian calendar.
      if (abs(ioptn) .ge. 3) then
!
!        Julian calendar values.
         yrref  =    1600
         jdref  = 2305518
!               = Julian Day reference value for the day that begins
!                 at noon on 1 March of the reference year "yrref".
         ndy400 = 400*365 + 100
         ndy100 = 100*365 +  25
!
!        Adjust for Gregorian calendar values.
         if (ioptn .gt. 0) then
            jdref  = jdref  - 10
            ndy400 = ndy400 -  3
            ndy100 = ndy100 -  1
         endif
      endif
!
!----------------------------------------------------------------
! OPTIONS -1 and +1:
! Given a calendar date (iday,month,iyear), compute the day number
! of the year (idayct), where 1 January is day number 1 and 31 December
! is day number 365 or 366, depending on whether it is a leap year.
      if (abs(ioptn) .eq. 1) then
!
!     Compute the day number during the year.
      if (month .le. 2) then
         idayct = iday + (month-1)*31
      else
         idayct = iday + int(30.6001 * (month+1)) - 63 + leap
      endif
!
!----------------------------------------------------------------
! OPTIONS -2 and +2:
! Given the day number of the year (idayct) and the year (iyear),
! compute the day of the month (iday) and the month (month).
      elseif (abs(ioptn) .eq. 2) then
!
      if (idayct .lt. 60+leap) then
         month  = (idayct-1)/31
         iday   = idayct - month*31
         month  = month + 1
      else
         ndays  = idayct - (60+leap)
!               = number of days past 1 March of the current year.
         jmonth = (10*(ndays+31))/306 + 3
!               = month counter, =4 for March, =5 for April, etc.
         iday   = (ndays+123) - int(30.6001*jmonth) 
         month  = jmonth - 1
      endif
!
!----------------------------------------------------------------
! OPTIONS -3 and +3:
! Given a calendar date (iday,month,iyear), compute the Julian Day
! number (idayct) that starts at noon.
      elseif (abs(ioptn) .eq. 3) then
!
!     Shift to a system where the year starts on 1 March, so January
!     and February belong to the preceding year.
!     Define jmonth=4 for March, =5 for April, ..., =15 for February.
      if (month .le. 2) then
        jyear  = jyear -  1
        jmonth = month + 13
      else
        jmonth = month +  1
      endif
!
!     Find the closest multiple of 400 years that is .le. jyear.
      yr400 = (jyear/400)*400
!           = multiple of 400 years at or less than jyear.
      if (jyear .lt. yr400) then
         yr400 = yr400 - 400
      endif
!
      n400yr = (yr400 - yrref)/400
!            = number of 400-year periods from yrref to yr400.
      nyrs   = jyear - yr400
!            = number of years from the beginning of yr400
!              to the beginning of jyear.
!
!     Compute the Julian Day number.
      idayct = iday + int(30.6001*jmonth) - 123 + 365*nyrs + nyrs/4 &
     &       + jdref + n400yr*ndy400
!
!     If we are using the Gregorian calendar, we must not count
!     every 100-th year as a leap year.  nyrs is less than 400 years,
!     so we do not need to consider the leap year that would occur if
!     nyrs were divisible by 400, i.e., we do not add nyrs/400.
      if (ioptn .gt. 0) then
         idayct = idayct - nyrs/100
      endif
!
!----------------------------------------------------------------
! OPTIONS -5, -4, +4, and +5:
! Given the Julian Day number (idayct) that starts at noon,
! compute the corresponding calendar date (iday,month,iyear)
! (abs(ioptn)=4) or day number during the year (abs(ioptn)=5).
      else
!
!     Create a new reference date which begins on the nearest
!     400-year cycle less than or equal to the Julian Day for 1 March
!     in the year in which the given Julian Day number (idayct) occurs.
      ndays  = idayct - jdref
      n400yr = ndays / ndy400
!            = integral number of 400-year periods separating
!              idayct and the reference date, jdref.
      jdref  = jdref + n400yr*ndy400
      if (jdref .gt. idayct) then
         n400yr = n400yr - 1
         jdref  = jdref  - ndy400
      endif
!
      ndays  = idayct - jdref
!            = number from the reference date to idayct.
!
      n100yr = min(ndays/ndy100, 3)
!            = number of complete 100-year periods
!              from the reference year to the current year.
!              The min() function is necessary to avoid n100yr=4
!              on 29 February of the last year in the 400-year cycle.
!
      ndays  = ndays - n100yr*ndy100
!            = remainder after removing an integral number of
!              100-year periods.
!
      n4yr   = ndays / 1461
!            = number of complete 4-year periods in the current century.
!              4 years consists of 4*365 + 1 = 1461 days.
!
      ndays  = ndays - n4yr*1461
!            = remainder after removing an integral number
!              of 4-year periods.
!
      n1yr   = min(ndays/365, 3)
!            = number of complete years since the last leap year.
!              The min() function is necessary to avoid n1yr=4
!              when the date is 29 February on a leap year,
!              in which case ndays=1460, and 1460/365 = 4.
!
      ndays  = ndays - 365*n1yr
!            = number of days so far in the current year,
!              where ndays=0 on 1 March.
!
      iyear  = n1yr + 4*n4yr + 100*n100yr + 400*n400yr + yrref 
!            = year, as counted in the standard way,
!              but relative to 1 March.
!
! At this point, we need to separate ioptn=abs(4), which seeks a
! calendar date, and ioptn=abs(5), which seeks the day number during
! the year.  First compute the calendar date if desired (abs(ioptn)=4).
      if (abs(ioptn) .eq. 4) then
         jmonth = (10*(ndays+31))/306 + 3
!               = offset month counter.  jmonth=4 for March, =13 for
!                 December, =14 for January, =15 for February.
         iday   = (ndays+123) - int(30.6001*jmonth)
!               = day of the month, starting with 1 on the first day
!                 of the month.
!
!        Now adjust for the fact that the year actually begins
!        on 1 January.
         if (jmonth .le. 13) then
            month = jmonth - 1
         else
            month = jmonth - 13
            iyear = iyear + 1
         endif
!
! This code handles abs(ioptn)=5, finding the day number during the year.
      else
!        ioptn=5 always returns month=1, which we set now.
         month = 1
!
!        We need to determine whether this is a leap year.
         leap = 0
         if ((jyear/4)*4 .eq. jyear) then
            leap = 1
         endif
         if ((ioptn .gt. 0)               .and. &
     &       ((jyear/100)*100 .eq. jyear) .and. &
     &       ((jyear/400)*400 .ne. jyear)      ) then
               leap = 0
         endif
!
!        Now find the day number "iday".
!        ndays is the number of days since the most recent 1 March,
!        so ndays=0 on 1 March.
         if (ndays .le.305) then
            iday  = ndays + 60 + leap
         else
            iday  = ndays - 305
            iyear = iyear + 1
         endif
      endif
!
!     Adjust the year if it is .le. 0, and hence BC (Before Christ).
      if (iyear .le. 0) then
         iyear = iyear - 1
      endif
!
! End the code for the last option, ioptn.
      endif
!
      return
      end


      integer function idaywk(jdayno)
!
! IDAYWK = compute the DAY of the WeeK given the Julian Day number,
!          version 1.0.
!
! Input variable
      integer  jdayno
! jdayno = Julian Day number starting at noon of the day in question.
!
! Output variable:
! idaywk = day of the week, where 0=Sunday, 1=Monday, ..., 6=Saturday.
!
!----------
! Compute the day of the week given the Julian Day number.
! You can find the Julian Day number given (day,month,year)
! using subroutine calndr.f.
! Example: For the first day of the Gregorian calendar,
! 15 October 1582, compute the Julian day number (option 3 of
! subroutine calndr) and compute the day of the week.
!     call calndr (3, 15, 10, 1582, jdayno) 
!     write(*,*) jdayno, idaywk(jdayno)
! The numbers printed should be 2299161 and 5,
! where 6 refers to Friday.
!
! Copyright (C) 1999 Jon Ahlquist.
! Issued under the second GNU General Public License.
! See www.gnu.org for details.
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
! If you find any errors, please notify:
! Jon Ahlquist <ahlquist@met.fsu.edu>   
! Dept of Meteorology
! Florida State University
! Tallahassee, FL 32306-4520
! 15 March 1999.
!
!-----
! Declare internal variable.
! jdSun is the Julian Day number starting at noon on any Sunday.
! I arbitrarily chose the first Sunday after Julian Day 1,
! which is Julian Day 6.
      integer  jdSun
      data     jdSun /6/
      idaywk = mod(jdayno-jdSun,7)
! If jdayno-jdSun < 0, then we are taking the modulus of a negative
! number. Fortran's built-in mod function returns a negative value
! when the argument is negative.  In that case, we adjust the result
! to a positive value.
      if (idaywk .lt. 0) idaywk = idaywk + 7
      return

      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

