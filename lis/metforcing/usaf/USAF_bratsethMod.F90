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
! MODULE: USAF_bratsethMod
!
! REVISION HISTORY:
! 22 Jun 2017  Initial version.........................Eric Kemp/SSAI/NASA
! 31 Aug 2018  Added IMERG support, plus some refactoring, code optimization,
!              and revised summaries of QC checks......Eric Kemp/SSAI/NASA
! 07 Sep 2018  Renamed to USAF_bratsethMod.F90.........Eric Kemp/SSAI/NASA
! 21 Feb 2020  Added support for 10-km GALWEM..........Eric Kemp/SSAI/NASA
! 05 Mar 2020  Added support for new GFS filename convention
!              ........................................Eric Kemp/SSAI/NASA
! 03 Jun 2020  Removed Box-Cox transform in precipitation analysis
!              ........................................Eric Kemp/SSAI/NASA
!
! DESCRIPTION:
!
! Source code for gridded precipitation and screen-level analyses using
! the Bratseth objective analysis scheme.  Includes routines for storing
! observations in data structures, partitioning precipitation into 3-hour
! increments, interpolating data, performing quality control, and running
! the analysis.  Some concepts are borrowed from the Canadian Precipitation
! Analysis and from research performed by ECMWF.
!
! USEFUL REFERENCES:
!
! Bratseth, A M, 1986:  Statistical interpolation by means of successive
!   corrections.  Tellus, 38A, 439-447.
! Cressie, N A C, 1993:  Statistics for Spatial Data.  Revised Edition, 
!   Wiley, New York, 928 pp.
! Daley, R, 1991:  Atmospheric Data Analysis.  Cambridge University Press,
!   Cambridge, UK, 457 pp.
! Fortin, V, G Roy, N Donaldson, and A Majidjiba, 2015:  Assimilation of
!   radar quantitative precipitation estimations in the Canadian
!   Precipitation Analysis (CaPA).  J Hydrol, 531, 296-307.
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

#include "LIS_misc.h"

! Macros for tracing - Requires ESMF 7_1_0+
#ifdef ESMF_TRACE
#define TRACE_ENTER(region) call ESMF_TraceRegionEnter(region)
#define TRACE_EXIT(region) call ESMF_TraceRegionExit(region)
#else
#define TRACE_ENTER(region)
#define TRACE_EXIT(region)
#endif

module USAF_bratsethMod

   ! Imports
#ifdef ESMF_TRACE
   use ESMF
#endif

   ! Defaults
   implicit none
   private

   ! Data structure for observed and background values at observation
   ! points.
   type USAF_obsData
      private
      integer :: nobs
      character*32, allocatable :: net(:)
      character*32, allocatable :: platform(:)
      real, allocatable :: obs(:) ! Observed variable
      real, allocatable :: lat(:) ! Latitude of observation (deg N)
      real, allocatable :: lon(:) ! Longitude of observation (deg E)
      real, allocatable :: sigmaOSqr(:) ! Error variance of observation
      real, allocatable :: oErrScaleLength(:) ! Obs error correlation length
      real, allocatable :: back(:) ! Background variable
      integer, allocatable :: qc(:) ! Quality control flag
   end type USAF_obsData
   public :: USAF_obsData

   ! Public methods
   public :: USAF_createObsData
   public :: USAF_destroyObsData
   public :: USAF_countGoodObs
   public :: USAF_assignObsData
   public :: USAF_multObsData
   public :: USAF_getBackNWP
   public :: USAF_split6hrGaugeObsData
   public :: USAF_split12hrGaugeObsData
   public :: USAF_getSSMIObsData
   public :: USAF_getGeoPrecipObsData
   public :: USAF_interpBackToTypeObsData
   public :: USAF_analyzePrecip
   public :: USAF_analyzeScreen
   public :: USAF_getCMORPHObsData
   public :: USAF_setBratsethPrecipStats
   public :: USAF_setBratsethScreenStats
   ! EMK New
   public :: USAF_filterObsData
   public :: USAF_waterQC
   public :: USAF_dupQC
   public :: USAF_snowQC
   public :: USAF_snowDepthQC
   public :: USAF_backQC
   public :: USAF_superstatQC

   ! A simple linked list type that can be used in a hash table.  Intended
   ! to store indices of arrays in the USAF_obsData type for efficient look-up.
   type hash_list
      integer :: obindex
      type(hash_list), pointer :: next
   end type hash_list

   ! Private constants
   real, parameter :: MISSING = -9999

   ! Quality control flags.  Should these be public?
   integer, parameter :: QC_UNKNOWN = 0
   integer, parameter :: QC_GOOD = 1
   integer, parameter :: QC_SUSPECT = 2
   integer, parameter :: QC_REJECT = 3

contains

   !---------------------------------------------------------------------------
   ! Constructor
   ! FIXME: Merge with initObsData
   subroutine USAF_createObsData(this,n,maxobs)

      ! Imports
      use AGRMET_forcingMod, only : agrmet_struc

      ! Defaults
      implicit none

      ! Arguments
      type(USAF_ObsData),intent(inout) :: this
      integer, intent(in) :: n
      integer, optional, intent(in) :: maxobs

!      TRACE_ENTER("bratseth_create")
      ! Allocate array members of ObsData type
      this%nobs = 0
      if (present(maxobs)) then
         allocate(this%net(maxobs))
         allocate(this%platform(maxobs))
         allocate(this%obs(maxobs))
         allocate(this%lat(maxobs))
         allocate(this%lon(maxobs))
         allocate(this%sigmaOSqr(maxobs))
         allocate(this%oErrScaleLength(maxobs))
         allocate(this%back(maxobs))
         allocate(this%qc(maxobs))
      else
         allocate(this%net(agrmet_struc(n)%max_pcpobs))
         allocate(this%platform(agrmet_struc(n)%max_pcpobs))
         allocate(this%obs(agrmet_struc(n)%max_pcpobs))
         allocate(this%lat(agrmet_struc(n)%max_pcpobs))
         allocate(this%lon(agrmet_struc(n)%max_pcpobs))
         allocate(this%sigmaOSqr(agrmet_struc(n)%max_pcpobs))
         allocate(this%oErrScaleLength(agrmet_struc(n)%max_pcpobs))
         allocate(this%back(agrmet_struc(n)%max_pcpobs))
         allocate(this%qc(agrmet_struc(n)%max_pcpobs))
      end if
      call initObsData(this)
!      TRACE_EXIT("bratseth_create")

   end subroutine USAF_createObsData

   !---------------------------------------------------------------------------
   ! Destructor
   subroutine USAF_destroyObsData(this)

      ! Defaults
      implicit none

      ! Arguments
      type(USAF_ObsData),intent(inout) :: this

!      TRACE_ENTER("bratseth_dstroy")
      ! Dellocate array members of ObsData type
      this%nobs = 0
      deallocate(this%net)
      deallocate(this%platform)
      deallocate(this%obs)
      deallocate(this%lat)
      deallocate(this%lon)
      deallocate(this%sigmaOSqr)
      deallocate(this%oErrScaleLength)
      deallocate(this%back)
      deallocate(this%qc)
!      TRACE_EXIT("bratseth_dstroy")

   end subroutine USAF_destroyObsData

   !---------------------------------------------------------------------------
   ! FIXME: Merge with USAF_createObsData
   subroutine initObsData(this)

      ! Defaults
      implicit none

      ! Arguments
      type(USAF_ObsData),intent(inout) :: this

      ! Initialize data members
      this%nobs = 0
      this%net(:) = "NULL"
      this%platform(:) = "NULL"
      this%obs(:) = MISSING
      this%lat(:) = MISSING
      this%lon(:) = MISSING
      this%sigmaOSqr(:) = MISSING
      this%oErrScaleLength(:) = MISSING
      this%back(:) = MISSING
      this%qc(:) = QC_UNKNOWN

   end subroutine InitObsData

   !---------------------------------------------------------------------------
   ! Loops through contents of ObsData structure and counts all obs with
   ! "good" quality control flag.
   function USAF_countGoodObs(this) result(goodObs)
      implicit none
      type(USAF_ObsData), intent(in) :: this
      integer :: goodObs
      integer :: i
!      TRACE_ENTER("bratseth_cntGood")
      goodObs = 0
      do i = 1, this%nobs
         if (this%qc(i) .ne. QC_REJECT) then
            goodObs = goodObs + 1
         end if
      end do ! i
!      TRACE_EXIT("bratseth_cntGood")
   end function USAF_countGoodObs

   !---------------------------------------------------------------------------
   ! Copies a single observation into a ObsData structure.  Value of
   ! background field at observation is optional (useful for adding
   ! "superobservations").
   subroutine USAF_assignObsData(this,net,platform,ob,lat,lon,sigmaOSqr, &
        oErrScaleLength,back)

      ! Imports
      use LIS_logmod, only : LIS_logunit

      ! Defaults
      implicit none

      ! Arguments
      type(USAF_ObsData),intent(inout) :: this
      character(len=32), intent(in) :: net
      character(len=32), intent(in) :: platform
      real, intent(in) :: ob
      real, intent(in) :: lat
      real, intent(in) :: lon
      real, intent(in) :: sigmaOSqr
      real, intent(in) :: oErrScaleLength
      real, optional, intent(in) :: back

      ! Local variables
      integer :: nobs

      ! Sanity check.  Since this is intended for an operational system,
      ! just print a warning and return if we see an array bounds problem.
      nobs = this%nobs
      if (nobs .eq. size(this%obs,1)) then
         write(LIS_logunit,*) &
              '[WARN], not enough memory for assigning observation data!'
         write(LIS_logunit,*) &
              '[WARN] Can only store ',nobs,' observations'
         return
      end if
!      TRACE_ENTER("bratseth_assign")

      ! Assign the value.
      nobs = nobs + 1
      this%net(nobs) = net
      this%platform(nobs) = platform
      this%obs(nobs) = ob
      this%lat(nobs) = lat
      ! Make sure -180 to 180 convention is used
      if (lon .gt. 180) then
         this%lon(nobs) = lon - 360.
      else
         this%lon(nobs) = lon
      end if
      this%sigmaOSqr(nobs) = sigmaOSqr
      this%oErrScaleLength(nobs) = oErrScaleLength
      this%qc(nobs) = QC_UNKNOWN
      this%nobs = nobs
      if (present(back)) then
         this%back(nobs) = back
      end if
!      TRACE_EXIT("bratseth_assign")

   end subroutine USAF_assignObsData

   !---------------------------------------------------------------------------
   ! Multiply the observed and background field values by a given value.
   ! Useful for unit conversions.
   subroutine USAF_multObsData(this,multiple)

      ! Defaults
      implicit none

      ! Arguments
      type(USAF_ObsData), intent(inout) :: this
      real, intent(in) :: multiple

      ! Local variables
      integer :: i

      TRACE_ENTER("bratseth_mltply")
      do i = 1, this%nobs
         this%obs(i) = multiple*this%obs(i)
         if (this%back(i) .eq. MISSING) cycle
         this%back(i) = multiple*this%back(i)
      end do ! i
      TRACE_EXIT("bratseth_mltply")
   end subroutine USAF_multObsData

   !---------------------------------------------------------------------------
   ! Disaggregate 6-hr precipitation accumulations into 3-hr periods using
   ! the background field.  New 3-hr values are stored in new ObsData
   ! structures.
   subroutine USAF_split6hrGaugeObsData(this,nest,imax,jmax,back4,pcap,p3,p6)

      ! Imports

      ! Defaults
      implicit none

      ! Arguments
      type(USAF_ObsData),intent(in) :: this
      integer,intent(in) :: nest
      integer,intent(in) :: imax
      integer,intent(in) :: jmax
      real, intent(in) :: back4(imax,jmax,4)
      real,intent(in) :: pcap
      type(USAF_ObsData),intent(out) :: p3
      type(USAF_ObsData),intent(out) :: p6

      ! Local variables 
      integer :: nobs6
      real :: tmp_back2(2)
      real :: tmp_obs2(2)
      real :: e6
      integer :: n,k
      real, allocatable :: backObsPts(:,:)
      integer, external :: get_fieldpos
      real, parameter :: FILL = MISSING

      TRACE_ENTER("bratseth_split06hr")

      ! See if we have any observations to process
      nobs6 = this%nobs
      if (nobs6 .eq. 0) then
          call USAF_createObsData(p3,nest,maxobs=1)
          call USAF_createObsData(p6,nest,maxobs=1)
         TRACE_EXIT("bratseth_split06hr")
         return
      end if
      !call USAF_createObsData(p3,nest,maxobs=nobs6)
      !call USAF_createObsData(p6,nest,maxobs=nobs6)
      call USAF_createObsData(p3,nest,maxobs=size(this%net))
      call USAF_createObsData(p6,nest,maxobs=size(this%net))

      ! Interpolate the background field to the observations for each
      ! time slice.
      allocate(backObsPts(nobs6,2))
      call interpBack2ObsPts(this,nest,imax,jmax,back4(:,:,1), &
           backObsPts(:,1))
      call interpBack2ObsPts(this,nest,imax,jmax,back4(:,:,2), &
           backObsPts(:,2))

      ! Divide 6-hr gauge reports into 3-hr increments based on background
      ! field.
      do n = 1, nobs6

         ! Skip for missing data
         if ( this%qc(n) .eq. QC_REJECT .or. &
              backObsPts(n,1) .eq. MISSING .or. &
              backObsPts(n,2) .eq. MISSING) then
            cycle
         end if

         tmp_back2(1:2) = backObsPts(n,1:2)
         tmp_obs2(1:2) = MISSING

         ! Get sum of two 3-hr background values
         e6 = 0.
         do k = 1,2
            e6 = e6 + tmp_back2(k)
         end do ! k

         ! If no 3-hr precip available from background, just divide
         ! the 6-hr gauge data equally.
         if (nint(e6) .eq. 0) then
            tmp_obs2(1:2) = 0.5 * this%obs(n)
         else
            ! Use the ratio of 3-hr background values to divide the gauge
            ! data
            do k = 1,2
               tmp_obs2(k) = (tmp_back2(k)/e6) * this%obs(n)
            end do
         end if

         ! Filter out bad obs
         do k = 1,2
            if (tmp_obs2(k) .gt. pcap) then
               tmp_obs2(k) = MISSING
            end if
         end do

         ! Assign 
         if (tmp_obs2(1) .ne. MISSING) then
            call USAF_assignObsData(p3, this%net(n), this%platform(n), &
                 tmp_obs2(1), this%lat(n), this%lon(n), this%sigmaOSqr(n), &
                 this%oErrScaleLength(n),back=tmp_back2(1))
         end if
         if (tmp_obs2(2) .ne. MISSING) then
            call USAF_assignObsData(p6, this%net(n), this%platform(n), &
                 tmp_obs2(2), this%lat(n), this%lon(n), this%sigmaOSqr(n), &
                 this%oErrScaleLength(n),back=tmp_back2(2))
         end if

      end do ! n
      TRACE_EXIT("bratseth_split06hr")
   end subroutine USAF_split6hrGaugeObsData

   !---------------------------------------------------------------------------
   ! Disaggregate 12-hr precipitation accumulations into 3-hr periods using
   ! the background field.  New 3-hr values are stored in new ObsData
   ! structures.
   subroutine USAF_split12hrGaugeObsData(this,nest,imax,jmax,back4,pcap,p3,p6,&
        p9,p12)

      ! Imports

      ! Defaults
      implicit none

      ! Arguments
      type(USAF_ObsData),intent(in) :: this
      integer,intent(in) :: nest
      integer,intent(in) :: imax
      integer,intent(in) :: jmax
      real, intent(in) :: back4(imax,jmax,4)
      real,intent(in) :: pcap
      type(USAF_ObsData),intent(out) :: p3
      type(USAF_ObsData),intent(out) :: p6
      type(USAF_ObsData),intent(out) :: p9
      type(USAF_ObsData),intent(out) :: p12

      ! Local variables
      integer :: nobs12
      real :: tmp_back4(4)
      real :: tmp_obs4(4)
      real :: e12
      integer :: n,k
      real, allocatable :: backObsPts(:,:)
      integer, external :: get_fieldpos
      real, parameter :: FILL = MISSING

      TRACE_ENTER("bratseth_split12hr")

      ! See if we have any observations to process
      nobs12 = this%nobs
      if (nobs12 .eq. 0) then
         call USAF_createObsData(p3,nest,maxobs=1)
         call USAF_createObsData(p6,nest,maxobs=1)
         call USAF_createObsData(p9,nest,maxobs=1)
         call USAF_createObsData(p12,nest,maxobs=1)
         TRACE_EXIT("bratseth_split12hr")
         return
      end if

      ! Initialize the 3, 6, 9, and 12 hour precipitation objects
      ! TODO:  Allow maxobs to be set via lis.config.
      !call USAF_createObsData(p3,nest,maxobs=nobs12)
      !call USAF_createObsData(p6,nest,maxobs=nobs12)
      !call USAF_createObsData(p9,nest,maxobs=nobs12)
      !call USAF_createObsData(p12,nest,maxobs=nobs12)
      call USAF_createObsData(p3,nest,maxobs=size(this%net))
      call USAF_createObsData(p6,nest,maxobs=size(this%net))
      call USAF_createObsData(p9,nest,maxobs=size(this%net))
      call USAF_createObsData(p12,nest,maxobs=size(this%net))

      ! Interpolate the background field to the observations for each
      ! time slice.
      allocate(backObsPts(nobs12,4))
      call interpBack2ObsPts(this,nest,imax,jmax,back4(:,:,1), &
           backObsPts(:,1))
      call interpBack2ObsPts(this,nest,imax,jmax,back4(:,:,2), &
           backObsPts(:,2))
      call interpBack2ObsPts(this,nest,imax,jmax,back4(:,:,3), &
           backObsPts(:,3))
      call interpBack2ObsPts(this,nest,imax,jmax,back4(:,:,4), &
           backObsPts(:,4))

      ! Divide 12-hr gauge reports into 3-hr increments based on background
      ! field.
      do n = 1, nobs12
         ! Skip for missing data
         if ( this%qc(n) .eq. QC_REJECT .or. &
              backObsPts(n,1) .eq. MISSING .or. &
              backObsPts(n,2) .eq. MISSING .or. &
              backObsPts(n,3) .eq. MISSING .or. &
              backObsPts(n,4) .eq. MISSING) then
            cycle
         end if

         tmp_back4(1:4) = backObsPts(n,1:4)
         tmp_obs4(1:4) = MISSING

         ! Get sum of four 3-hr background values
         e12 = 0.
         do k = 1,4
            e12 = e12 + tmp_back4(k)
         end do ! k

         ! If no 3-hr precip available from background, just divide
         ! the 12-hr gauge data equally.
         if (nint(e12) .eq. 0) then
            tmp_obs4(1:4) = 0.25 * this%obs(n)
         else
            ! Use the ratio of 3-hr background values to divide the gauge
            ! data
            do k = 1,4
               tmp_obs4(k) = (tmp_back4(k)/e12) * this%obs(n)
            end do
         end if

         ! Filter out bad obs
         do k = 1,4
            if (tmp_obs4(k) .gt. pcap) then
               tmp_obs4(k) = MISSING
            end if
         end do

         ! Assign
         if (tmp_obs4(1) .ne. MISSING) then
            call USAF_assignObsData(p3, this%net(n), this%platform(n), &
                 tmp_obs4(1), this%lat(n), this%lon(n), this%sigmaOSqr(n), &
                 this%oErrScaleLength(n),back=tmp_back4(1))
         end if
         if (tmp_obs4(2) .ne. MISSING) then
            call USAF_assignObsData(p6, this%net(n), this%platform(n), &
                 tmp_obs4(2), this%lat(n), this%lon(n), this%sigmaOSqr(n), &
                 this%oErrScaleLength(n),back=tmp_back4(2))
         end if
         if (tmp_obs4(3) .ne. MISSING) then
            call USAF_assignObsData(p9, this%net(n), this%platform(n), &
                 tmp_obs4(3), this%lat(n), this%lon(n), this%sigmaOSqr(n), &
                 this%oErrScaleLength(n),back=tmp_back4(3))
         end if
         if (tmp_obs4(4) .ne. MISSING) then
            call USAF_assignObsData(p12, this%net(n), this%platform(n), &
                 tmp_obs4(4), this%lat(n), this%lon(n), this%sigmaOSqr(n), &
                 this%oErrScaleLength(n),back=tmp_back4(4))
         end if

      end do ! n
      TRACE_EXIT("bratseth_split12hr")
   end subroutine USAF_split12hrGaugeObsData

   !---------------------------------------------------------------------------
   ! Copy gridded SSMI precipitation retrievals into ObsData structure.
   ! Map projection logic borrowed from AGRMET_forcingMod.F90.
   ! ONLY POLAR STEREOGRAPHIC 8TH MESH AND 16TH MESH GRIDS ARE SUPPORTED.
   subroutine USAF_addSSMIObsData(this,imax,jmax,ra_tmp,nest)

      ! Imports
      use AGRMET_forcingMod, only : agrmet_struc
      use LIS_logMod, only: LIS_logunit, LIS_endrun

      ! Defaults
      implicit none

      ! Arguments
      type(USAF_ObsData), intent(inout) :: this
      integer,intent(in) :: imax
      integer,intent(in) :: jmax
      real,intent(in) :: ra_tmp(2,imax,jmax)
      integer, intent(in) :: nest

      ! Local variables
      integer :: ihemi
      integer :: npts
      real, allocatable :: xpts(:), ypts(:), rlat(:), rlon(:)
      real :: sigmaOSqr, ob, xi1, xj1, oErrScaleLength
      real :: xpnmcaf, ypnmcaf, orient, xmesh, xmeshl
      character(len=32) :: net
      character(len=32) :: platform
      integer :: icount
      integer :: i,j
      integer :: count_good_ssmi
      real, parameter :: FILL = -9999.0

      external :: polarToLatLon

      net = "SSMI"
      platform = "SSMI"
      sigmaOSqr = agrmet_struc(nest)%bratseth_precip_ssmi_sigma_o_sqr
      oErrScaleLength = &
           agrmet_struc(nest)%bratseth_precip_ssmi_err_scale_length

      ! Allocate memory for map projection calculations
      npts = imax*jmax
      allocate(xpts(npts))
      allocate(ypts(npts))
      allocate(rlat(npts))
      allocate(rlon(npts))

      ! Need to assemble SSMI map projection based on hemisphere and number
      ! of grid points.  Taken from AGRMET_forcingMod.F90
      ! SSM/I imax is either 512, 1024, or 1440.
      if (imax .eq. 512) then
         xmeshl = 47.625
         xpnmcaf = 257
         ypnmcaf = 257
      else if (imax .eq. 1024) then
         xmeshl = 47.625/2
         xpnmcaf = 513
         ypnmcaf = 513
      else if (imax .eq. 1440) then
         continue
      else
         write(LIS_logunit,*)'[ERR] Invalid imax dimension for SSM/I!'
         write(LIS_logunit,*)'[ERR] Received ', imax
         write(LIS_logunit,*)'[ERR] Only support 512, 1024, or 1400'
         flush(LIS_logunit)
         call LIS_endrun()
      end if

      ! Handle the polar stereographic case first
      if (imax .eq. 512 .or. imax .eq. 1024) then
         count_good_ssmi = 0
         do ihemi = 1, 2

            ! Assemble gridDesc array
            if(ihemi.eq.1) then
               xmesh = xmeshl
               orient = 100.0
            else
               xmesh = -xmeshl
               orient = 280.0
            endif

            ! Get lat/lon of SSM/I data.
            icount = 1
            do j = 1,jmax
               xj1 = real(j) - ypnmcaf
               do i = 1,imax
                  xi1 = float(i) - xpnmcaf
                  call polarToLatLon(xi1,xj1,xmesh,orient,&
                       rlat(icount),rlon(icount))
                  icount = icount + 1
               end do ! i
            end do ! j

            ! At this point, we have the lat/lon for SSM/I in this hemisphere
            ! and the LIS global map projection.
            icount = 0
            do j = 1, jmax
               do i = 1, imax

                  icount = icount + 1

                  ! Screen out water points
                  if (agrmet_struc(nest)%land(i,j,ihemi) .eq. 0) cycle

                  ! Latitude bounds checks for hemispheric grids
                  if (ihemi .eq. 1) then
                     if (rlat(icount) .lt. 0) cycle
                  else if (ihemi .eq. 2) then
                     if (rlat(icount) .gt. 0) cycle
                  end if

                  if (ra_tmp(ihemi,i,j) .lt. 0) cycle
                  ob = ra_tmp(ihemi,i,j)
                  count_good_ssmi = count_good_ssmi + 1
                  call USAF_assignObsData(this,net,platform,ob, rlat(icount), &
                       rlon(icount),sigmaOSqr,oErrScaleLength)
               end do ! i
            end do ! j

         end do ! ihemi

         write(LIS_logunit,*) &
              '[INFO] Found ',count_good_ssmi,' SSMI obs'
      end if  ! Polar stereographic case

      ! FIXME Handle Lat/Lon case
      if (imax .eq. 1440) then

         write(LIS_logunit,*)'[ERR] Lat/lon SSM/I data not supported yet!'
         write(LIS_logunit,*)'[ERR] Modify USAF_addSSMIObsData and recompile!'
         flush(LIS_logunit)
         call LIS_endrun()

      end if ! Lat/Lon case

      ! Clean up
      deallocate(xpts)
      deallocate(ypts)
      deallocate(rlat)
      deallocate(rlon)

   end subroutine USAF_addSSMIObsData

   !---------------------------------------------------------------------------
   ! Interpolate background precipitation fields to LIS grid and collect
   ! in global domain array.  This is salvaged from older AGRMET code
   ! for legacy Barnes analysis.  Supports GFS and GALWEM.
   subroutine USAF_getBackNWP(nest,back4,pcp_src, use_twelve, j6hr, findex)

      ! Imports
      use AGRMET_forcingMod, only: agrmet_struc
      use LIS_coreMod, only: LIS_rc, LIS_masterproc
      use LIS_logMod, only: LIS_abort, LIS_endrun, LIS_getNextUnitNumber, &
           LIS_releaseUnitNumber, LIS_logunit, LIS_alert
      use LIS_mpiMod
      use LIS_timeMgrMod, only: LIS_julhr_date

      ! Defaults
      implicit none

      ! Arguments
      integer, intent(in) :: nest
      real, intent(inout) :: back4(LIS_rc%gnc(nest), LIS_rc%gnr(nest), 4)
      character(len=6), intent(inout) :: pcp_src(4)
      logical, intent(in) :: use_twelve
      integer, intent(in) :: j6hr
      integer, intent(in) :: findex

      ! Local variables
      real :: fg_data_glb(LIS_rc%gnc(nest), LIS_rc%gnr(nest))
      integer :: k, t
      integer :: j3hr
      integer :: fc_hr
      character(len=6) :: src
      character(len=255) :: message(20)
      integer :: rc,ierr
      character(len=10) :: yyyymmddhh
      integer :: c, r

      external :: AGRMET_julhr_date10

      TRACE_ENTER("bratseth_getBackNWP")
      rc = 0

      ! For (0, 6Z], k = 1,2.  For (6,12Z], k = 3,4
      if (use_twelve) then
         k = 2
         do t=3,4
            back4(:,:,t) = LIS_rc%udef ! EMK
         end do
      else
         k = 0
         do t=1,2
            back4(:,:,t) = LIS_rc%udef ! EMK
         end do
      end if

      fg_data_glb = LIS_rc%udef

      do j3hr = j6hr+3, j6hr+6, 3
         k = k + 1
         ierr = 0

         call AGRMET_julhr_date10(j3hr, yyyymmddhh)

         ! EMK...Fetch bias ratio if requested
         if (agrmet_struc(nest)%back_bias_corr .eq. 1) then
            call USAF_pcpBackBiasRatio_s2s(nest, yyyymmddhh)
         else if (agrmet_struc(nest)%back_bias_corr .eq. 2) then
            call USAF_pcpBackBiasRatio_nrt(nest, yyyymmddhh)
         end if

         write(LIS_logunit,*) &
              '[INFO] Searching for NWP precipitation valid ', yyyymmddhh

         ! Get GALWEM, if requested
         if (agrmet_struc(nest)%galwemprecswch .eq. 1) then
            src = "GALWEM"
            if ((k .eq. 1) .or. (k .eq. 3)) then
               fc_hr=3
               call fldbld_precip_nwp(nest,findex,j6hr,src,fc_hr, &
                    fg_data_glb,rc)
               ierr = rc
            else if ((k .eq. 2) .or. (k .eq. 4)) then
               fc_hr=6
               call fldbld_precip_nwp(nest,findex,j6hr,src,fc_hr,&
                    fg_data_glb,rc)
               ierr = rc
            end if

            ! Use GALWEM Bratseth settings if we have the data
            if (rc .eq. 0) then
               pcp_src(k) = 'GALWEM'

               ! Apply bias correction to background field
               if (agrmet_struc(nest)%back_bias_corr .eq. 2) then
                  write(LIS_logunit,*) &
           '[INFO] Applying IMERG-based bias correction to GALWEM precip'
                  do r = 1, LIS_rc%gnr(nest)
                     do c = 1, LIS_rc%gnc(nest)
                        fg_data_glb(c,r) = &
                             fg_data_glb(c,r) * &
                          agrmet_struc(nest)%galwem_nrt_bias_ratio(c,r)
                     end do
                  end do
               end if

            end if
         endif

         ! Get GFS, if requested, or if we have no GALWEM data
         if (agrmet_struc(nest)%gfsprecswch .eq. 1 .or. &
              (ierr .ne. 0)) then

            if (ierr .ne. 0) then
               write(LIS_logunit,*) &
                    '[WARN] Falling back on GFS precipitation... '
            end if

            ierr = 0
            src = "GFS"
            if ((k .eq. 1) .or. (k .eq. 3)) then
               fc_hr=3
               call fldbld_precip_nwp(nest,findex,j6hr,src,fc_hr, &
                    fg_data_glb,rc)
               ierr = rc
            else if ((k .eq. 2) .or. (k .eq. 4)) then
               fc_hr=6
               call fldbld_precip_nwp(nest,findex,j6hr,src,fc_hr, &
                    fg_data_glb,rc)
               ierr = rc
            end if

            ! Apply bias correction to background field
            if (agrmet_struc(nest)%back_bias_corr .eq. 1 .and. &
                 rc .eq. 0) then
               write(LIS_logunit,*) &
             '[INFO] Applying GALWEM-based bias correction to GFS precip'
               do r = 1, LIS_rc%gnr(nest)
                  do c = 1, LIS_rc%gnc(nest)
                     fg_data_glb(c,r) = &
                          fg_data_glb(c,r) * &
                          agrmet_struc(nest)%pcp_back_bias_ratio(c,r)
                  end do
               end do

            else if (agrmet_struc(nest)%back_bias_corr .eq. 2 .and. &
                 rc .eq. 0) then
               ! New IMERG option for NRT
               write(LIS_logunit,*) &
             '[INFO] Applying IMERG-based bias correction to GFS precip'
               do r = 1, LIS_rc%gnr(nest)
                  do c = 1, LIS_rc%gnc(nest)
                     fg_data_glb(c,r) = &
                          fg_data_glb(c,r) * &
                          agrmet_struc(nest)%gfs_nrt_bias_ratio(c,r)
                  end do
               end do
            end if

            ! Use GFS Bratseth settings if we have the data
            if (rc .eq. 0) then
               pcp_src(k) = "GFS   "
            end if
         endif

         ! Handle missing precip
         if (ierr .ne. 0) then
            call AGRMET_julhr_date10(j3hr, yyyymmddhh)
            write(LIS_logunit,*) &
                 '[ERR] No NWP background precipitation found for ',yyyymmddhh
            write(LIS_logunit,*) '[ERR] ABORTING!'
            flush(LIS_logunit)
            message(:) = ''
            message(1) = '[ERR] Program:  LIS'
            message(2) = '  Routine:  USAF_getBackNWP.'
            message(3) = '  GRIB precipitation data not available for '//&
                 yyyymmddhh

#if (defined SPMD)
            call MPI_Barrier(LIS_MPI_COMM, ierr)
            call handle_mpi_error(ierr, &
                 'MPI_Barrier call in USAF_getBackNWP')
#endif
            if(LIS_masterproc) then
               call LIS_alert( 'LIS.USAF_getBackNWP          ', 1, &
                    message )
               call LIS_abort( message)
            else
               call sleep(10) ! Make sure LIS_masterproc finishes LIS_abort
               call LIS_endrun()
            endif
         end if

         back4(:,:,k) = fg_data_glb(:,:)

      end do ! j3hr
      TRACE_EXIT("bratseth_getBackNWP")

   end subroutine USAF_getBackNWP

   !---------------------------------------------------------------------------
   ! Get gridded SSMI rainfall retrievals, and copy to appropriate 3-hr
   ! ObsData structure.  Calls USAF_addSSMIObsData under the hood.  Borrows
   ! logic from earlier AGRMET code for legacy Barnes analysis.
   subroutine USAF_getSSMIObsData(nest,j6hr,use_twelve,precip3,precip6, &
        precip9,precip12,pcp_src)

      ! Imports
      use AGRMET_forcingMod, only : agrmet_struc
      use LIS_coreMod, only: LIS_masterproc
      use LIS_fileIOMod, only: LIS_putget
      use LIS_logMod, only: LIS_logunit, LIS_alert
      use LIS_timeMgrMod, only : LIS_julhr_date

      ! Defaults
      implicit none

      ! Arguments
      integer, intent(in) :: nest
      integer, intent(in) :: j6hr
      logical, intent(in) :: use_twelve
      type(USAF_obsData), intent(inout) :: precip3
      type(USAF_obsData), intent(inout) :: precip6
      type(USAF_obsData), intent(inout) :: precip9
      type(USAF_obsData), intent(inout) :: precip12
      character(len=6), intent(in) :: pcp_src(4)

      ! Local variables
      real, allocatable :: ra(:,:,:)
      integer :: k
      integer :: j3hr
      integer :: imax,jmax
      integer :: alert_number
      integer :: first,last,hemi
      integer :: yr,mo,da,hr
      character(len=100) :: ifil
      logical :: exists
      character(len=255) :: message(20)
      character(len=30) :: routine_name
      logical :: use_zeros
      integer :: local_global_or_hemi
      integer :: local_diff_grid
      integer :: ii,jj,kk,count_good_ssmi ! EMK TEST

      external :: agrmet_ssmiprec_filename

      TRACE_ENTER("bratseth_getSSMI")
      routine_name = 'USAF_getSSMIObsData'

      alert_number = 0
      if (use_twelve) then
         k = 2
      else
         k = 0
      end if

      ! Loop through time levels
      do j3hr = j6hr+3, j6hr+6, 3
         k = k + 1

         ! Set Bratseth error statistics based on source of background field.
         call USAF_setBratsethPrecipStats(pcp_src(k),nest)

         if (agrmet_struc(nest)%raswch .eq. 1) then
            local_global_or_hemi = agrmet_struc(nest)%global_or_hemi
            if (agrmet_struc(nest)%imaxsmi .eq. 1440) then
               local_global_or_hemi=1
            end if
            local_diff_grid = agrmet_struc(nest)%diff_grid
            if ( agrmet_struc(nest)%imaxsmi .ne. &
                 agrmet_struc(nest)%imax) then
               local_diff_grid=1
            end if

            ! Select appropriate SSMI grid dimensions
            if (local_diff_grid .eq. 1) then
               imax = agrmet_struc(nest)%imax3
               jmax = agrmet_struc(nest)%jmax3
            else
               imax = agrmet_struc(nest)%imax
               jmax = agrmet_struc(nest)%jmax
            end if

            allocate(ra(2,imax,jmax))
!            ra = LIS_rc%udef
            ra = MISSING
            use_zeros = .false.
            if (agrmet_struc(nest)%razero .eq. 1) then
               use_zeros = .true.
            end if

            ! Get the data
            if (local_global_or_hemi.eq.0) then
!               first=agrmet_struc(nest)%shemi
!               last=agrmet_struc(nest)%nhemi
               first = 1
               last = 2
            else
               first=1
               last=1
            endif
            do hemi=first,last
                call LIS_julhr_date( j3hr, yr,mo,da,hr)
                call agrmet_ssmiprec_filename(ifil,&
                     agrmet_struc(nest)%agrmetdir,&
                     agrmet_struc(nest)%ssmidir,&
                     agrmet_struc(nest)%use_timestamp,&
                     hemi,yr,mo,da,hr,nest,imax,jmax)

                inquire(file = trim(ifil), exist = exists)
                if (.not. exists) then
                   write(LIS_logunit,*) ' '
                   write(LIS_logunit,*) &
                        '[WARN] precip/smiedr:  error opening file'
                   write(LIS_logunit,*)'[WARN]  SSMI data file ', trim(ifil), &
                        ' does not exist.'
                   write(LIS_logunit,*)'[WARN]  SSMI estimates will not be used ',&
                        'in precip analysis.'
                   write(LIS_logunit,*) ' '
                   message   =' '
                   message(1)='[WARN] program:  LIS'
                   message(2)='  routine:  AGRMET_smiest'
                   message(3)='  SSMI data file '//trim(ifil)//&
                        ' does not exist.'
                   message(4)='  SSMI estimates will not be used in '//&
                        'precip analysis.'
                   alert_number = alert_number + 1
                   if(LIS_masterproc) &
                        call lis_alert('precip              ', alert_number, &
                        message )

                   ra(hemi,:,:) = MISSING
                else
                   write(LIS_logunit,*) '[INFO] READING ',trim(ifil)
                   call LIS_putget( ra(hemi,:,:), 'r', ifil, routine_name, &
                        imax, jmax)
                end if ! .not. exists
            end do ! hemi

            ! Honor option to reset SSMI zero precip values to missing
            if (.not. use_zeros) then
               write(LIS_logunit,*)'[INFO] SSMI ZEROS NOT USED'
               where ( .not. ra(:,:,:) > 0.0 .and. &
                    .not. ra(:,:,:) < 0)
                  ra(:,:,:) = MISSING
               end where
            end if
            where ( ra(:,:,:) .eq. 9999)
               ra(:,:,:) = MISSING
            end where

            count_good_ssmi = 0
            do jj = 1,jmax
               do ii = 1,imax
                  do kk = 1,2
                     if (ra(kk,ii,jj) .gt. 0) then
                        count_good_ssmi = count_good_ssmi + 1
                     end if
                  end do
               end do
            end do
            write(LIS_logunit,*)'[INFO] Found ',count_good_ssmi, ' SSMI obs'

            ! Now append the SSMI values to the appropriate 3-hr precip
            ! structure.  We will interpolate background values later.
            if (k .eq. 1) then
               call USAF_addSSMIObsData(precip3,imax,jmax,ra,nest)
            else if (k .eq. 2) then
               call USAF_addSSMIObsData(precip6,imax,jmax,ra,nest)
            else if (k .eq. 3) then
               call USAF_addSSMIObsData(precip9,imax,jmax,ra,nest)
            else
               call USAF_addSSMIObsData(precip12,imax,jmax,ra,nest)
            end if

            ! Clean up
            deallocate(ra)

         end if ! raswch
      end do ! j3hr
      TRACE_EXIT("bratseth_getSSMI")
   end subroutine USAF_getSSMIObsData

   !---------------------------------------------------------------------------
   ! Get gridded GEOPRECIP rainfall retrievals, and copy to appropriate 3-hr
   ! ObsData structure.  Calls USAF_addSSMIObsData under the hood.  Borrows
   ! logic from earlier AGRMET code for legacy Barnes analysis.
   subroutine USAF_getGeoPrecipObsData(nest,j6hr,use_twelve,precip3,precip6,&
        precip9,precip12,pcp_src)

      ! Imports
      use AGRMET_forcingMod, only : agrmet_struc
      use LIS_coreMod, only     : LIS_masterproc
      use LIS_fileIOMod, only: LIS_putget
      use LIS_logMod, only: LIS_logunit, LIS_endrun, LIS_alert
      use LIS_timeMgrMod, only : LIS_julhr_date

      ! Defaults
      implicit none

      ! Arguments
      integer, intent(in) :: nest
      integer, intent(in) :: j6hr
      logical, intent(in) :: use_twelve
      type(USAF_obsData), intent(inout) :: precip3
      type(USAF_obsData), intent(inout) :: precip6
      type(USAF_obsData), intent(inout) :: precip9
      type(USAF_obsData), intent(inout) :: precip12
      character(len=6),intent(in) :: pcp_src(4)

      ! Local variables
      real, allocatable :: geoprc(:,:)
      real, allocatable :: geopr(:,:)
      integer :: local_global_or_hemi
      integer :: local_diff_grid
      integer :: j3hr
      integer :: mesh
      integer :: imax,jmax
      integer :: first, last, hemi
      logical :: gdgeornk, exists
      integer :: yr,mo,da,hr
      character(len=100) :: ifil
      integer :: k
      integer :: alert_number
      character(len=255) :: message(20)
      character(len=30) :: routine_name
      real, allocatable :: xpts(:), ypts(:), rlat(:), rlon(:)
      real :: sigmaOSqr, ob, xi1, xj1, oErrScaleLength
      real :: xpnmcaf, ypnmcaf, orient, xmesh, xmeshl
      character(len=32) :: net
      character(len=32) :: platform
      integer :: count_good_geo_precip, icount
      integer :: npts
      integer :: i,j,jj
      logical, external :: is_geo_corrupted
      real, allocatable :: gest_temp(:,:,:)

      external :: agrmet_geoprec_filename
      external :: polarToLatLon

      TRACE_ENTER("bratseth_getGeopPrcp")
      net = "GEOPRECIP"
      platform = "GEOPRECIP"
      sigmaOSqr = agrmet_struc(nest)%bratseth_precip_geoprecip_sigma_o_sqr
      oErrScaleLength = &
           agrmet_struc(nest)%bratseth_precip_geoprecip_err_scale_length

      routine_name = "USAF_getGeoPrecipObsData"
      alert_number = 0

      if (use_twelve) then
         k = 2
      else
         k = 0
      end if

      ! Loop through time levels
      do j3hr = j6hr+3, j6hr+6, 3
         k = k + 1

         ! Set Bratseth error statistics based on source of background field.
         call USAF_setBratsethPrecipStats(pcp_src(k),nest)

         if ( (agrmet_struc(nest)%geoswch .eq. 1) .or. &
              (agrmet_struc(nest)%geoswch .eq. 2)) then
            local_global_or_hemi = agrmet_struc(nest)%global_or_hemi
            if (agrmet_struc(nest)%imaxgp == 3600) then
               local_global_or_hemi = 1
            end if
            local_diff_grid = agrmet_struc(nest)%diff_grid
            if (agrmet_struc(nest)%imaxgp .ne. agrmet_struc(nest)%imax) then
               local_diff_grid = 1
            end if

            ! Select mesh number and array sizes
            mesh = 1
            imax = agrmet_struc(nest)%imax
            jmax = agrmet_struc(nest)%jmax
            if (local_diff_grid .eq. 1) then
               imax = agrmet_struc(nest)%imax2
               jmax = agrmet_struc(nest)%jmax2
               if (agrmet_struc(nest)%imaxgp .eq. 4096) then
                  if (agrmet_struc(nest)%imax .eq. 512) then
                     mesh = 8
                  else if (agrmet_struc(nest)%imax .eq. 1024) then
                     mesh = 4
                  else
                     write(LIS_logunit,*) &
                          '[ERR] Invalid dimension for GEO_PRECIP data!'
                     write(LIS_logunit,*)'[ERR] Read ', agrmet_struc(nest)%imax
                     call LIS_endrun()
                  end if
               end if
            end if

            ! Allocate memory for map projection calculations
            npts = imax*jmax
            allocate(xpts(npts))
            allocate(ypts(npts))
            allocate(rlat(npts))
            allocate(rlon(npts))

            ! Assemble GEO_PRECIP map projection.  Taken from &
            ! AGRMET_forcingMod.F90
            if (imax .eq. 512) then
               xmeshl = 47.625
               xpnmcaf = 257
               ypnmcaf = 257
            else if (imax .eq. 1024) then
               xmeshl = 47.625/2
               xpnmcaf = 513
               ypnmcaf = 513
            else if (imax .eq. 4096) then
               xmeshl = 47.625/8
               xpnmcaf = 2049
               ypnmcaf = 2049
            else
               write(LIS_logunit,*) &
                    '[ERR] Invalid imax dimension for GEO_PRECIP!'
               write(LIS_logunit,*)'[ERR] Received ',imax
               write(LIS_logunit,*)'[ERR] Only support 512, 1024, or 4096!'
               flush(LIS_logunit)
               call LIS_endrun()
            end if

            if (local_global_or_hemi .eq. 0) then
!               first = agrmet_struc(nest)%shemi
!               last  = agrmet_Struc(nest)%nhemi
               first = 1
               last = 2
            else
               first = 1
               last  = 1
            endif

            do hemi = first,last

               gdgeornk = .false.
               call LIS_julhr_date( j3hr, yr,mo,da,hr)
               call agrmet_geoprec_filename(ifil,agrmet_struc(nest)%agrmetdir,&
                    agrmet_struc(nest)%geodir, &
                    agrmet_struc(nest)%use_timestamp,&
                    hemi,yr,mo,da,hr,nest,imax,jmax)

               inquire( file = trim(ifil), exist = exists)
               if ( .not. exists ) then
                  write(LIS_logunit,*) &
                       '[WARN] USAF_getGeoPrecipObsData:  error opening file ', &
                       trim(ifil)
                  write(LIS_logunit,*) '[WARN]  file does not exist'
                  write(LIS_logunit,*) &
                       '[WARN]  geo precip estimate will not be performed'
                  message = ' '
                  message(1) = 'program:  LIS'
                  message(2) = '  routine:  USAF_getGeoPrecipObsData'
                  message(3) = '  error opening file '//trim(ifil)
                  message(4) = '  file does not exist'
                  message(5) = '  geo precip estimate will not be performed'
                  alert_number = alert_number + 1
                  if(LIS_masterproc) then
                     call LIS_alert( 'LIS.USAF_getGeoPrecipObsData', &
                          alert_number, message )
                  endif
                  TRACE_EXIT("bratseth_getGeopPrcp")
                  return
               endif
               write(LIS_logunit,*) '[INFO] READING ', trim(ifil)

               allocate(geoprc(imax,jmax))

               call LIS_putget( geoprc, 'r', ifil, routine_name, &
                    imax, jmax )

               ! Need to shift 10th degree GEO_PRECIP by 180 degrees; the data
               ! start at 0.05, not at -179.95
               if (imax .eq. 3600) then
                  allocate(geopr(imax,jmax))
                  geopr=cshift(geoprc,1800,dim=1)
                  geoprc=geopr
                  deallocate(geopr)
               end if

               ! Check for anomalous geoprecip files
               if (is_geo_corrupted(geoprc, imax, jmax, mo, hemi)) then
                  write(LIS_logunit,*) &
                       '[WARN] USAF_getGeoPrecipObsData:  data corrupted - ', &
                       trim(ifil)
                  write(LIS_logunit,*) &
                       '[WARN]  geo precip estimate will not be performed'
                  message = ' '
                  message(1) = 'program:  LIS'
                  message(2) = '  routine:  USAF_getGeoPrecipObsData'
                  message(3) = '  data corrupted in file '//trim(ifil)
                  message(4) = '  '
                  message(5) = '  geo precip estimate will be assigned missing'
                  alert_number = alert_number + 1
                  if(LIS_masterproc) then
                     call lis_alert( 'LIS.USAF_getGeoPrecipObsData', &
                          alert_number, message )
                  endif
                  geoprc = MISSING
               end if

               ! Loop through field and reset bad points to MISSING
               do j = 1, jmax
                  do i = 1, imax
                     if (geoprc(i,j) .lt. 0) then
                        geoprc(i,j)  = MISSING
                     end if
                  end do
               end do
               if ( local_diff_grid .eq. 1 .and. &
                    agrmet_struc(nest)%imaxgp .ne. 4096) then
                  do j = 1, jmax
                     do i = 1, imax
                        if (agrmet_struc(nest)%land2(i,j,hemi) .le. 0) then
                           geoprc(i,j) = MISSING
                        else
                           if (geoprc(i,j) .ge. 9990.0) then
                              geoprc(i,j) = MISSING
                           end if
                        end if
                     end do ! i
                  end do ! j
               else
                  do j = 1, jmax
                     do i = 1, imax
                        if (agrmet_struc(nest)%land( &
                             int(i/mesh),int(j/mesh),hemi) .le. 0) then
                           geoprc(i,j) = MISSING
                        else
                           if (geoprc(i,j) .ge. 9990.0) then
                              geoprc(i,j) = MISSING
                           end if
                        end if
                     end do ! i
                  end do ! j
               end if

               ! Flip the data if necessary
               if (local_global_or_hemi .eq. 1) then
                  allocate(gest_temp(1,imax,jmax))
                  do j = 1, jmax
                     jj = jmax-(j-1)
                     do i = 1, imax
                        gest_temp(1,i,j) = geoprc(i,jj)
                     end do ! i
                  end do ! j
                  geoprc(:,:) = gest_temp(1,:,:)
                  deallocate(gest_temp)
               end if

               ! Handle the polar stereographic case first
               if ( imax .eq. 512 .or. imax .eq. 1024 .or. &
                    imax .eq. 4096) then
                  if (hemi .eq. 1) then
                     xmesh = xmeshl
                     orient = 100.0
                  else
                     xmesh = -xmeshl
                     orient = 280.0
                  end if

                  ! Get lat/lon of GEO_PRECIP data.
                  icount = 1
                  do j = 1, jmax
                     xj1 = real(j) - ypnmcaf
                     do i = 1,imax
                        xi1 = float(i) - xpnmcaf
                        call polarToLatLon(xi1,xj1,xmesh,orient, &
                             rlat(icount),rlon(icount))
                        icount = icount + 1
                     end do ! i
                  end do ! j

                  ! At this point, we have the lat/lon for GEO_PRECIP in this
                  ! hemisphere, plus the LIS global map projection.
                  icount = 0
                  count_good_geo_precip = 0
                  do j = 1, jmax
                     do i = 1, imax

                        icount = icount + 1

                        ! Skip GEO_PRECIP if over water
                        if (agrmet_struc(nest)%land(i,j,hemi) .eq. 0) cycle

                        ! Lat/Lon bounds
                        if (rlat(icount) .gt.  50.) cycle
                        if (rlat(icount) .lt. -50.) cycle

                        if (geoprc(i,j) .eq. MISSING) cycle

                        ob = geoprc(i,j)
                        count_good_geo_precip = count_good_geo_precip + 1
                        if (k .eq. 1) then
                           call USAF_assignObsData(precip3,net,platform,ob, &
                                rlat(icount),rlon(icount),sigmaOSqr, &
                                oErrScaleLength)
                        else if (k .eq. 2) then
                           call USAF_assignObsData(precip6,net,platform,ob, &
                                rlat(icount),rlon(icount),sigmaOSqr, &
                                oErrScaleLength)
                        else if (k .eq. 3) then
                           call USAF_assignObsData(precip9,net,platform,ob, &
                                rlat(icount),rlon(icount),sigmaOSqr, &
                                oErrScaleLength)
                        else if (k .eq. 4) then
                           call USAF_assignObsData(precip12,net,platform,ob, &
                                rlat(icount),rlon(icount),sigmaOSqr, &
                                oErrScaleLength)
                        end if
                     end do ! i
                  end do ! j

                  write(LIS_logunit,*)'[INFO] count_good_geo_precip = ', &
                       count_good_geo_precip

               end if ! Polar stereographic cases

               ! Handle lat/lon case
               if (imax .eq. 3600) then

                  write(LIS_logunit,*)&
                       '[ERR] Lat/lon GEO_PRECIP data not supported yet!'
                  write(LIS_logunit,*)&
                       '[ERR] Modify USAF_addGeoPrecipObsData and recompile!'
                  flush(LIS_logunit)
                  call LIS_endrun()

               end if

               ! Clean up
               deallocate(geoprc)

            end do ! hemi

            deallocate(rlat)
            deallocate(rlon)
            deallocate(xpts)
            deallocate(ypts)

         end if ! geoswch
      end do ! j3hr
      TRACE_EXIT("bratseth_getGeopPrcp")

   end subroutine USAF_getGeoPrecipObsData

   !---------------------------------------------------------------------------
   ! Interpolates background field to selected observation points, and
   ! appends the interpolated values to the ObsData structure.  Useful for
   ! working with superobservations.
   subroutine USAF_interpBackToTypeObsData(this,nest,imax,jmax,back,type)

      ! Imports

      ! Defaults
      implicit none

      ! Arguments
      type(USAF_ObsData),intent(inout) :: this
      integer,intent(in) :: nest
      integer,intent(in) :: imax
      integer,intent(in) :: jmax
      real, intent(in) :: back(imax,jmax)
      character(len=32),intent(in) :: type

      ! Local variables
      integer :: nobs
      integer :: n
      real, allocatable :: backObsPts(:)
      integer, external :: get_fieldpos
      real, parameter :: FILL = MISSING

      ! See if we have any observations to process
      nobs = this%nobs
      if (nobs .eq. 0) return

!      TRACE_ENTER("bratseth_interp2Typ")

      ! Interpolate the background field to the observations for the selected
      ! time slice.
      allocate(backObsPts(nobs))
      call interpBack2ObsPts(this,nest,imax,jmax,back, &
           backObsPts)

      ! Loop through the observations and insert the background values to
      ! the observation type.
      do n = 1, nobs

         if (.not. trim(type) .eq. "ALL") then
            if ( trim(this%net(n)) .ne. trim(type) .and. &
                 trim(this%platform(n)) .ne. trim(type)) cycle
         end if

         this%back(n) = backObsPts(n)

      end do ! n

      ! Cleanup
      deallocate(backObsPts)
!      TRACE_EXIT("bratseth_interp2Typ")

   end subroutine USAF_interpBackToTypeObsData

   !---------------------------------------------------------------------------
   ! Surface precipitation analysis using a NWP background field and
   ! irregularly positioned observations.
   !
   ! (1) Currently supported input data are rain gauges, CMORPH estimates,
   !     SSMI retrievals, GEOPRECIP retrievals, and IMERG estimates.
   ! (2) The observations are subjected to quality control tests (all
   !     before calling analyzePrecip subroutine).
   ! (3) All observations (including superobservations) that pass quality
   !     control are merged into a single ObsData structure.
   ! (4) The inverse data density is calculated for each observation (see
   !     Bratseth 1986).
   ! (5) The Bratseth analysis is generated via iteration at each observation
   !     point.  This successive correction scheme converges to Optimal
   !     Interpolation without direct matrix inversion, providing significant
   !     computational savings and allowing a pseudo-global analysis (the
   !     radius of influence is specified by the background error covariance
   !     rather than arbitrary limits on the number of observations).
   !     The observed, background, and analysis values are also collected in
   !     a OBA structure for later output to file.
   ! (6) The Bratseth analysis is interpolated to the LIS grid in a single
   !     pass, using an algorithm similar to those of Daley (1991)
   !     Pedder (1993), and Kalnay (2003).
   !
   ! NOTES:
   ! (1) This implementation uses a Gaussian function to model error
   !     covariances.  Semivariogram analyses suggest the inverse exponential
   !     function has a better fit to actual statistics, but the Gaussian
   !     function has a much shorter radius of influence that greatly speeds
   !     up the analysis.
   ! (2) The IMERG product offers the future potential for useful frozen
   !     precipitation estimates as well as estimates over snow.  But the
   !     fidelity of these estimates are heavily dependent on the instrument
   !     used (e.g., GPM high frequency imager channels; sounder channels).  At
   !     present the code is conservative and rejects likely snowy
   !     observations, but this could be reconsidered in the future.
   ! (3) Other bias corrections and quality control tests for the input data
   !     are pending.
   !---------------------------------------------------------------------------

   subroutine USAF_analyzePrecip(precipAll,nest,back,hourindex,mrgp,precipOBA)

      ! Imports
      use AGRMET_forcingMod, only:  agrmet_struc
      use LIS_coreMod, only: LIS_rc, LIS_ews_halo_ind, LIS_ewe_halo_ind, &
           LIS_nss_halo_ind, LIS_nse_halo_ind, LIS_localPet
      use USAF_OBAMod, only: OBA, createOBA,assignOBA

      ! Defaults
      implicit none

      ! Arguments
      type(USAF_ObsData),intent(inout) :: precipAll
      integer,intent(in) :: nest
      real, intent(in) :: back(LIS_rc%gnc(nest), LIS_rc%gnr(nest))
      integer, intent(in) :: hourindex
      real, intent(inout) :: mrgp(LIS_rc%lnc(nest),LIS_rc%lnr(nest))
      type(OBA), intent(out) :: precipOBA

      ! Local variables
      integer :: nobs
      real, allocatable :: invDataDensities(:)
      real, allocatable :: sumObsEstimates(:)
      integer :: npasses
      real :: convergeThresh
      real :: sigmaBSqr
      integer :: good_obs
      integer :: j

      TRACE_ENTER("bratseth_analyzePrcp")

      ! Initialize merged field with the background first guess.  This will be
      ! changed below as needed.
      mrgp(:,:) = back(LIS_ews_halo_ind(nest,LIS_localPet+1): &
           LIS_ewe_halo_ind(nest,LIS_localPet+1), &
           LIS_nss_halo_ind(nest,LIS_localPet+1): &
           LIS_nse_halo_ind(nest,LIS_localPet+1))

      ! If no observations, skip the Bratseth analysis
      if ( precipAll%nobs .eq. 0) then
         if (agrmet_struc(nest)%oba_switch .eq. 1 .or. &
              agrmet_struc(nest)%oba_switch .eq. 2) then
            precipOBA = createOBA(nest,maxobs=0)
         end if
         call USAF_destroyObsData(precipAll)
         call zeroTrace(LIS_rc%lnc(nest),LIS_rc%lnr(nest),mrgp)
         TRACE_EXIT("bratseth_analyzePrcp")
         return
      end if
      nobs = precipAll%nobs

      ! EMK...Option 2 just captures O and B info, and skips the
      ! analysis.
      if (agrmet_struc(nest)%oba_switch .eq. 2) then
         good_obs = USAF_countGoodObs(precipAll)
         precipOBA = createOBA(nest, maxobs=good_obs)
         do j = 1, precipAll%nobs
            if (precipAll%qc(j) .eq. QC_REJECT) cycle
            call assignOBA(precipOBA, &
                 precipAll%net(j), precipAll%platform(j), &
                 precipAll%lat(j), precipAll%lon(j), &
                 precipAll%obs(j), precipAll%back(j), &
                 A=0.)
         end do ! j
         call USAF_destroyObsData(precipAll)
         TRACE_EXIT("bratseth_analyzePrcp")
         return
      end if

      ! Calculate (inverse) data density around each observation.
      sigmaBSqr = agrmet_struc(nest)%bratseth_precip_back_sigma_b_sqr
      call calc_invDataDensities(precipAll,sigmaBSqr,nest, &
           agrmet_struc(nest)%bratseth_precip_max_dist, &
           agrmet_struc(nest)%bratseth_precip_back_err_scale_length, &
           is_gauge, &
           invDataDensities)

      ! Run Bratseth analysis at observation points, and collect the sum of
      ! the corrections at each observation point (in sumObsEstimates), along
      ! with the required number of iterations (npasses).  Also return
      ! OBA information for output.
      convergeThresh = 0.01
      call calc_obsAnalysis(precipAll,sigmaBSqr,nobs,invDataDensities,nest,&
           agrmet_struc(nest)%bratseth_precip_max_dist, &
           agrmet_struc(nest)%bratseth_precip_back_err_scale_length, &
           convergeThresh, is_gauge, sumObsEstimates, &
           npasses, precipOBA)

      ! Calculate analysis at grid points.
      call calc_gridAnalysis(precipAll,nest,sigmaBSqr,nobs,invDataDensities,&
        sumObsEstimates,npasses,back, &
        agrmet_struc(nest)%bratseth_precip_max_dist, &
        agrmet_struc(nest)%bratseth_precip_back_err_scale_length, &
        mrgp)

      ! Clean up
      deallocate(invDataDensities)
      deallocate(sumObsEstimates)
      call USAF_destroyObsData(precipAll)

      ! Clobber spurious negative or near-zero values.
      call reset_negative_values(nest,mrgp)
      call zeroTrace(LIS_rc%lnc(nest),LIS_rc%lnr(nest),mrgp)
      TRACE_EXIT("bratseth_analyzePrcp")

   end subroutine USAF_analyzePrecip

   ! **Private routines**
   !---------------------------------------------------------------------------
   ! Low-level routine for interpolating gridded background field to
   ! observations.  Output data are in an array, and are *not* appended
   ! to the ObsData structure by this routine.
   subroutine interpBack2ObsPts(this,nest,imax,jmax,back, &
        backObsPts)

      ! Imports
      use LIS_coreMod, only: LIS_rc

      ! Defaults
      implicit none

      ! Arguments
      type(USAF_ObsData),intent(in) :: this
      integer,intent(in) :: nest
      integer,intent(in) :: imax
      integer,intent(in) :: jmax
      real, intent(in) :: back(imax,jmax)
      real, intent(out) :: backObsPts(this%nobs)

      ! Local variables
      integer :: nobs
      integer, allocatable :: n11(:), n12(:), n21(:), n22(:)
      real, allocatable :: w11(:), w12(:), w21(:), w22(:)
      logical*1 :: lb(imax*jmax)
      logical*1, allocatable :: lo(:)
      real, allocatable :: xpts(:), ypts(:)
      integer :: mi,mo,nv
      integer :: n
      integer :: i1,i2,j1,j2
      real ::  xi, xf, yi, yf
      integer :: iret
      real :: glbGridDesc(50)
      integer, external :: get_fieldpos
      real, parameter :: FILL = MISSING

      external :: compute_grid_coord
      external :: bilinear_interp

      ! See if we have any observations available
      nobs = this%nobs
      if (nobs .eq. 0) return

      ! Allocate internal arrays
      allocate(n11(nobs))
      allocate(n12(nobs))
      allocate(n21(nobs))
      allocate(n22(nobs))
      allocate(w11(nobs))
      allocate(w12(nobs))
      allocate(w21(nobs))
      allocate(w22(nobs))
      allocate(xpts(nobs))
      allocate(ypts(nobs))
      allocate(lo(nobs))

      ! Create a modified gridDesc that is for the entire LIS domain, not just
      ! the local process region.
      glbGridDesc(:) = LIS_rc%gridDesc(nest,:)
      glbGridDesc(2) = glbGridDesc(32) ! gnc
      glbGridDesc(3) = glbGridDesc(33) ! gnr
      glbGridDesc(4) = glbGridDesc(34) ! lat(1,1)
      glbGridDesc(5) = glbGridDesc(35) ! lon(1,1)
      glbGridDesc(7) = glbGridDesc(37) ! lat(gnc,gnr)
      glbGridDesc(8) = glbGridDesc(38) ! lon(gnc,gnr)

      ! Recalculate grid coordinates of LIS tile
      mo = nobs
      call compute_grid_coord(glbGridDesc,mo,FILL,xpts,ypts, &
           this%lon,this%lat,nv)

      ! Calculate corners and weights for each observation point.
      do n = 1,mo
         xi = xpts(n)
         yi = ypts(n)
         w11(n) = 0.
         w21(n) = 0.
         w12(n) = 0.
         w22(n) = 0.
         if (xi.ne.FILL .and. yi.ne.FILL) then
            i1=xi
            i2=i1+1
            j1=yi
            j2=j1+1
            xf=xi-i1
            yf=yi-j1
            n11(n)=get_fieldpos(i1,j1,glbGridDesc)
            n21(n)=get_fieldpos(i2,j1,glbGridDesc)
            n12(n)=get_fieldpos(i1,j2,glbGridDesc)
            n22(n)=get_fieldpos(i2,j2,glbGridDesc)
            if(min(n11(n),n21(n),n12(n),n22(n)).gt.0) then
               w11(n)=(1-xf)*(1-yf)
               w21(n)=xf*(1-yf)
               w12(n)=(1-xf)*yf
               w22(n)=xf*yf
            else
               n11(n)=0
               n21(n)=0
               n12(n)=0
               n22(n)=0
            endif
         else
            n11(n)=0
            n21(n)=0
            n12(n)=0
            n22(n)=0
         endif
      end do ! n

      ! Now interpolate
      mi = imax*jmax
      lb = .true.
      call bilinear_interp(glbGridDesc,lb,&
              back,lo,backObsPts, mi,mo, &
              this%lat, this%lon, w11,w12,w21,w22,n11,n12,n21,n22,&
              MISSING, iret)

      ! Clean up
      deallocate(n11)
      deallocate(n12)
      deallocate(n21)
      deallocate(n22)
      deallocate(w11)
      deallocate(w12)
      deallocate(w21)
      deallocate(w22)
      deallocate(xpts)
      deallocate(ypts)
      deallocate(lo)

   end subroutine interpBack2ObsPts

   !---------------------------------------------------------------------------
   ! Calculate (inverse) data density around each observation.  Part of
   ! Bratseth scheme.  See Bratseth (1986) or Sashegyi et al (1993).
   ! Note that this implementation accounts for possible correlated
   ! observation errors (e.g., for satellite data).
   !
   ! This implementation is parallelized by distributing work among the
   ! different LIS ESMF PETs.  To improve performance, a 2D hash table is
   ! constructed to group observations by LIS grid box; this allows the code
   ! to avoid ob-to-ob comparisons that are obviously too far away to be
   ! correlated.
   !
   ! NOTE:  Requires LIS to be run in lat-lon projection!
   subroutine calc_invDataDensities(this,sigmaBSqr,nest,max_dist, &
        backErrScaleLength,isUncorrObType,invDataDensities,silent)

      ! Imports
      use LIS_coreMod, only: LIS_localPet, LIS_rc
      use LIS_logMod, only : LIS_logunit, LIS_endrun, LIS_endrun
      use LIS_mpiMod

      ! Defaults
      implicit none

      ! Arguments
      type(USAF_ObsData), intent(inout) :: this
      real, intent(in) :: sigmaBSqr
      integer, intent(in) :: nest
      real, intent(in) :: max_dist
      real, intent(in) :: backErrScaleLength
      logical, external :: isUncorrObType
      real, allocatable, intent(out) :: invDataDensities(:)
      logical, intent(in), optional :: silent

      ! Local variables
      type(hash_list), allocatable, target :: hash2d(:,:)
      real, allocatable :: dataDensities_pet(:)
      integer :: nobs
      integer :: r,c,i,j,iob,job
      integer :: lr,ur,lc1,lc2,rc1,rc2
      real :: dist
      real :: b, num, denom
      integer, allocatable :: iobs_neighbors_vector(:), jobs_cr_vector(:)
      integer :: nobs_neighbors, nobs_cr
      double precision :: t1, t2
      integer :: ierr
      integer :: pet, pet_incr
      integer :: imax,jmax
      logical :: verbose

      verbose = .true.
      if (present(silent)) then
         if (silent) verbose = .false.
      end if

      nobs = this%nobs
      if (nobs .eq. 0) return

      if (verbose) then
         write(LIS_logunit,*) &
              '[INFO] Calculating local data densities for ',nobs,' obs...'
      end if

#if (defined SPMD)
      call MPI_Barrier(LIS_MPI_COMM, ierr)
      call handle_mpi_error(ierr, &
           'MPI_Barrier call in calc_invDataDensities')
      t1 = MPI_Wtime()
#endif

      ! Here we create a 2d hash table storing the index values of each ob
      ! in linked lists for each LIS grid box.  This can help us screen
      ! out obviously unnecessary ob comparisons later.
      call build_hash2d(this,nest,imax,jmax,hash2d)

      allocate(dataDensities_pet(nobs))
      dataDensities_pet(:) = 0

      ! Now we need to loop through the 2D hash, estimate the influence
      ! bounds for considering obs in neighboring grid boxes, and calculate the
      ! data densities.
      ! NOTE:  We're assuming LIS is run with a Lat-Lon grid.
      pet = -1
      pet_incr = 1
      do r = 1, LIS_rc%gnr(nest)
         do c = 1, LIS_rc%gnc(nest)

            ! Make sure obs are actually in this box.
            if (hash2d(c,r)%obindex .eq. MISSING) cycle

            ! See which PET is responsible for this grid box.
            call update_pet(pet,pet_incr)
            if (pet .ne. LIS_localPet) cycle

            ! Find neighbors with positive correlation in background field.
            call find_gridpt_neighbors(c,r,nest,max_dist, &
                 lr,ur,lc1,rc1,lc2,rc2)

            ! Get list of obs in current grid box
            call get_hash2d(LIS_rc%gnc(nest),LIS_rc%gnr(nest),hash2d,c,r,&
                 nobs_cr,jobs_cr_vector)
            if (nobs_cr .eq. 0) cycle

            ! Get list of obs from neighboring grid boxes
            call get_neighbor_obs(LIS_rc%gnc(nest),LIS_rc%gnr(nest),hash2d,&
                 lr,ur,lc1,rc1,lc2,rc2,&
                 nest, nobs_neighbors, iobs_neighbors_vector)

            ! For each ob in the current grid box, calculate data density
            ! contributions from neighbors.
            do j = 1, nobs_cr
               job = jobs_cr_vector(j)
               do i = 1, nobs_neighbors
                  iob = iobs_neighbors_vector(i)

                  if (iob .eq. job) then
                     dist = 0
                  else
                     dist = &
                          great_circle_distance(this%lat(iob), &
                          this%lon(iob), this%lat(job), this%lon(job))
                  end if
                  if (dist .gt. max_dist) cycle

                  b = backErrCov(sigmaBSqr,dist,backErrScaleLength)
                  num = b
                  if (iob .eq. job) then
                     num = num + this%sigmaOSqr(job)
                  else if (trim(this%net(iob)) .eq. trim(this%net(job))) then
                     ! Satellite observations have correlated errors.
                     if (.not. isUncorrObType(this%net(job))) then
                        if (.not. this%oErrScaleLength(job) > 0 .and. &
                             .not. this%oErrScaleLength(job) < 0) then
                           write(LIS_logunit,*) &
                                '[ERR]: job, network, oErrScaleLength: ', &
                                job, trim(this%net(job)), &
                                this%oErrScaleLength(job)
                           flush(LIS_logunit)
                        end if

                        num = num + &
                             obsErrCov(this%sigmaOSqr(job), &
                                       this%oErrScaleLength(job), &
                                       dist)
                     end if
                  end if

                  denom = &
                       (sigmaBSqr + this%sigmaOSqr(iob)) * &
                       (sigmaBSqr + this%sigmaOSqr(job))
                  denom = sqrt(denom)

                  dataDensities_pet(job) = &
                       dataDensities_pet(job) + (num/denom)

               end do ! i
            end do ! j

            ! Clean up
            deallocate(iobs_neighbors_vector)
            deallocate(jobs_cr_vector)

         end do ! c
      end do ! r

      ! Clean up
      call destroy_hash2d(LIS_rc%gnc(nest),LIS_rc%gnr(nest),hash2d)
      deallocate(hash2d)

      ! Collect the results
#if (defined SPMD)
      allocate(invDataDensities(nobs))
      invDataDensities(:) = 0
      call MPI_Barrier(LIS_MPI_COMM, ierr)
      call handle_mpi_error(ierr, &
           'MPI_Barrier call in calc_invDataDensities')
      call MPI_ALLREDUCE(dataDensities_pet,invDataDensities,nobs,MPI_REAL, &
           MPI_SUM, LIS_MPI_COMM, ierr)
      call handle_mpi_error(ierr, &
           'MPI_ALLREDUCE call in calc_invDataDensities')
#endif

      ! Clean up
      deallocate(dataDensities_pet)

      ! Finish data density calculations, including inversion.  This is
      ! probably fast enough to not warrant parallelization.
      do j = 1, nobs
         ! Skip bad data
         if ( this%qc(j) .eq. QC_REJECT) cycle
         invDataDensities(j) = &
              invDataDensities(j)*(sigmaBSqr + this%sigmaOSqr(j))
         invDataDensities(j) = 1. / invDataDensities(j)
      end do ! j

#if (defined SPMD)
      call MPI_Barrier(LIS_MPI_COMM, ierr)
      call handle_mpi_error(ierr, &
           'MPI_Barrier call in calc_invDataDensities')
      t2 = MPI_Wtime()
      if (verbose) then
         write(LIS_logunit,*) &
              '[INFO] Elapsed time calculating data densities is ',t2 - t1, &
              ' seconds'
      end if
#endif

   end subroutine calc_invDataDensities

   !---------------------------------------------------------------------------
   ! Perform Bratseth analysis at observation points.  Multiple iterations
   ! are made until convergence is reached.  Along the way, the observation
   ! estimates from each iteration are summed for later interpolation to the
   ! grid points.  Note that the *analysis* is also run at the observation
   ! points because in practice the analysis converges before the iterative
   ! "observation estimates" do (see Sashegyi et al 1993).
   !
   ! This implementation is parallelized by farming work out to the LIS ESMF
   ! PETs.  To improve performance, a 2D hash table is constructed to group
   ! observations by LIS grid box; this allows the code to avoid ob-to-ob
   ! comparisons that are obviously too far away to be correlated.
   !
   ! The observed, background, and analysis values at the observation
   ! points are also collected in an OBA structure for post-processing.
   subroutine calc_obsAnalysis(this,sigmaBSqr,nobs,invDataDensities,nest,&
        max_dist,backErrScaleLength, convergeThresh, isUncorrObType, &
        sumObsEstimates, npasses, varOBA, &
        skip, silent)

      ! Imports
      use AGRMET_forcingMod, only: agrmet_struc
      use LIS_coreMod, only : LIS_localPet, LIS_rc
      use LIS_logMod, only : LIS_logunit, LIS_endrun
      use LIS_mpiMod
      use USAF_OBAMod, only: OBA, createOBA, assignOBA

      ! Defaults
      implicit none

      ! Arguments
      type(USAF_ObsData), intent(inout) :: this
      real, intent(in) :: sigmaBSqr
      integer, intent(in) :: nobs
      real, intent(in) :: invDataDensities(nobs)
      integer, intent(in) :: nest
      real, intent(in) :: max_dist
      real, intent(in) :: backErrScaleLength
      real, intent(in) :: convergeThresh
      logical, external :: isUncorrObType
      real, allocatable, intent(out) :: sumObsEstimates(:)
      integer, intent(out) :: npasses
      type(OBA), intent(inout) :: varOBA
      logical, intent(in), optional :: skip(this%nobs)
      logical, intent(in), optional :: silent

      ! Local variables
      real, allocatable :: pprev_est(:)
      real, allocatable :: pprev_ana(:)
      real, allocatable :: pnew_ana(:)
      real, allocatable :: pnew_ana_pet(:)
      real, allocatable :: pnew_est(:)
      real, allocatable :: pnew_est_pet(:)
      real, allocatable :: sumObsEstimates_pet(:)
      type(hash_list), allocatable, target :: hash2d(:,:)
      integer :: cmax,rmax
      integer :: pet, pet_incr
      real :: dist, b, weight
      logical :: done
      integer :: imaxabsdiff
      real :: maxabsdiff, y_prev, y_new, normdev
      integer :: icount, num_high_dev
      integer :: c,r,i,j,iob,job,ii
      integer :: ierr
      double precision :: t0,t1, t2
      logical :: verbose
      real :: mad_est, mad_ana, diff
      real :: y_obs, y_est, y_ana
      integer :: lr,ur,lc1,rc1,lc2,rc2
      integer :: nobs_cr, nobs_neighbors
      integer, allocatable :: jobs_cr_vector(:), iobs_neighbors_vector(:)
      real :: O, A
      integer :: good_obs

      verbose = .true.
      if (present(silent)) then
         if (silent) verbose = .false.
      end if

      if (verbose) then
         write(LIS_logunit,*) &
              '[INFO] Running analysis at observation points...'
      end if

      ! Sanity checks
      if (nobs .eq. 0) return

      if (nobs .ne. this%nobs) then
         write(LIS_logunit,*) &
              '[ERR] Array size mismatch in calc_obsAnalysis!'
         write(LIS_logunit,*) &
              '[ERR] nobs, this%nobs = ',nobs, this%nobs
         call LIS_endrun()
      end if

      ! EMK TEST
      !do ii = 1, nobs
      !   write(LIS_logunit,*)'EMK: ii, invDataDensity: ', ii, &
      !        invDataDensities(ii)
      !end do

      ! Here we create a 2d hash table storing the index values of each ob
      ! in linked lists for each LIS grid box.  This can help us screen
      ! out obviously unnecessary ob comparisons later.
      call build_hash2d(this,nest,cmax,rmax,hash2d)

      ! Perform analysis at observation points.  See Bratseth (1986) or
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
      allocate(pnew_ana(nobs))
      allocate(pnew_ana_pet(nobs))
      allocate(pnew_est(nobs))
      allocate(pnew_est_pet(nobs))
      allocate(sumObsEstimates(nobs))
      allocate(sumObsEstimates_pet(nobs))
      allocate(pprev_ana(nobs))
      allocate(pprev_est(nobs))

      pprev_est(:) = this%back(:) ! First guess
      pprev_ana(:) = this%back(:) ! First guess
      sumObsEstimates(:) = 0
      sumObsEstimates_pet(:) = 0
      npasses = 0

#if (defined SPMD)
      call MPI_Barrier(LIS_MPI_COMM, ierr)
      call handle_mpi_error(ierr, &
           'MPI_Barrier call in calc_obsAnalysis')
      t0 = MPI_Wtime()
#endif

      do ! Iterate until convergence

#if (defined SPMD)
         call MPI_Barrier(LIS_MPI_COMM, ierr)
         call handle_mpi_error(ierr, &
           'MPI_Barrier call in calc_obsAnalysis')
         t1 = MPI_Wtime()
#endif
         pnew_est_pet(:) = 0
         pnew_ana_pet(:) = 0

         pet = -1
         pet_incr = 1

         ! Now we need to loop through the 2D hash, estimate the influence
         ! bounds for considering obs in neighboring grid boxes, and
         ! calculate analysis updates.
         do r = 1, LIS_rc%gnr(nest)
            do c = 1, LIS_rc%gnc(nest)

               ! Make sure obs are actually in this box.
               if (hash2d(c,r)%obindex .eq. MISSING) cycle

               ! See which PET is responsible for this grid box.
               call update_pet(pet,pet_incr)
               if (pet .ne. LIS_localPet) cycle

               ! Find neighbors with positive correlation in background field.
               call find_gridpt_neighbors(c,r,nest,max_dist, &
                    lr,ur,lc1,rc1,lc2,rc2)

               ! Get list of obs in current grid box
               call get_hash2d(LIS_rc%gnc(nest),LIS_rc%gnr(nest),hash2d,c,r,&
                    nobs_cr,jobs_cr_vector)
               if (nobs_cr .eq. 0) cycle

               ! Get list of obs from neighboring grid boxes
               call get_neighbor_obs(LIS_rc%gnc(nest),LIS_rc%gnr(nest),hash2d,&
                    lr,ur,lc1,rc1,lc2,rc2,&
                    nest, nobs_neighbors, iobs_neighbors_vector)

               ! For each ob in the current grid box, calculate analysis update
               ! contributions from neighbors.
               do j = 1, nobs_cr
                  job = jobs_cr_vector(j)

                  if (this%qc(job) .eq. QC_REJECT) then
                     sumObsEstimates_pet(j) = 0
                     cycle
                  endif

                  do i = 1, nobs_neighbors
                     iob = iobs_neighbors_vector(i)

                     if (this%qc(iob) .eq. QC_REJECT) then
                        sumObsEstimates_pet(iob) = 0
                        cycle
                     endif

                     if (iob .eq. job) then
                        dist = 0
                     else
                        dist = &
                             great_circle_distance(this%lat(iob), &
                             this%lon(iob), this%lat(job), this%lon(job))
                     end if
                     if (dist .gt. max_dist) cycle

                     b = backErrCov(sigmaBSqr,dist, &
                          backErrScaleLength)

                     ! First, update the observation estimate
                     weight = b
                     if (iob .eq. job) then
                        weight = weight + this%sigmaOSqr(iob)
                     else if (trim(this%net(iob)) .eq. &
                              trim(this%net(job))) then
                        ! Satellite data have horizontal error correlations
                        if (.not. isUncorrObType(this%net(job))) then
                           weight = weight + &
                                obsErrCov(this%sigmaOSqr(job), &
                                this%oErrScaleLength(job), &
                                dist)
                        end if
                     end if
                     weight = weight * invDataDensities(iob)
                     pnew_est_pet(job) = pnew_est_pet(job) + &
                          (weight*(this%obs(iob) - pprev_est(iob)))

                     ! Second, update the analysis at the observation point.
                     weight = b * invDataDensities(iob)
                     pnew_ana_pet(job) = pnew_ana_pet(job) + &
                          (weight*(this%obs(iob) - pprev_est(iob)))

                  end do ! i
               end do ! j
               if (allocated(iobs_neighbors_vector)) &
                    deallocate(iobs_neighbors_vector)
               if (allocated(jobs_cr_vector)) &
                    deallocate(jobs_cr_vector)
            end do ! c
         end do ! r

#if (defined SPMD)
         ! Share pnew_est and pnew_ana across all processors.
         pnew_est(:) = 0
         pnew_ana(:) = 0
         call MPI_Barrier(LIS_MPI_COMM, ierr)
         call handle_mpi_error(ierr, &
              'MPI_Barrier call in calc_obsAnalysis')
         call MPI_ALLREDUCE(pnew_est_pet, pnew_est, nobs, MPI_REAL, &
              MPI_SUM, LIS_MPI_COMM, ierr)
         call handle_mpi_error(ierr, &
              'MPI_ALLREDUCE call in calc_obsAnalysis')
         call MPI_Barrier(LIS_MPI_COMM, ierr)
         call handle_mpi_error(ierr, &
              'MPI_Barrier call in calc_obsAnalysis')
         call MPI_ALLREDUCE(pnew_ana_pet, pnew_ana, nobs, MPI_REAL, &
              MPI_SUM, LIS_MPI_COMM, ierr)
         call handle_mpi_error(ierr, &
              'MPI_ALLREDUCE call in calc_obsAnalysis')
#endif

         ! Finish analysis and observation estimates for this iteration
         do j = 1, nobs
            pnew_est(j) = pprev_est(j) + pnew_est(j)
            pnew_ana(j) = pprev_ana(j) + pnew_ana(j)
         end do ! r

         ! Update sum of observation estimates
         pet = -1
         pet_incr = 1
         do j = 1, nobs
            ! See which MPI process is responsible for this row.
            call update_pet(pet,pet_incr)
            if (pet .ne. LIS_localPet) cycle
            sumObsEstimates_pet(j) = sumObsEstimates_pet(j) + pprev_est(j)
         end do ! j
         npasses = npasses + 1

#if (defined SPMD)
         ! Share sumObsEstimates across all processors.
         sumObsEstimates(:) = 0
         call MPI_Barrier(LIS_MPI_COMM, ierr)
         call handle_mpi_error(ierr, &
              'MPI_Barrier call in calc_obsAnalysis')
         call MPI_ALLREDUCE(sumObsEstimates_pet, sumObsEstimates, nobs, &
              MPI_REAL, MPI_SUM, LIS_MPI_COMM, ierr)
         call handle_mpi_error(ierr, &
              'MPI_ALLREDUCE call in calc_obsAnalysis')
#else
         do j = 1, nobs
            sumObsEstimates(j) = sum(sumObsEstimates_pet)
         end do
#endif

         ! See if analysis has converged by comparing current and prior
         ! analysis values.  In practice, the analysis values will converge
         ! before the observation estimates do (unless observations are
         ! assumed perfect).  Then, overwrite prior values.
         done = .true.
         maxabsdiff = 0
         imaxabsdiff = 0
         normdev = 0
         num_high_dev = 0

         if (verbose) then
            mad_ana = 0
            mad_est = 0
            icount = 0
         end if

         do j = 1, nobs

            if (this%qc(j) .eq. QC_REJECT) cycle

            ! Check for convergence
            y_prev = pprev_ana(j)
            y_new = pnew_ana(j)
            if (abs(y_prev - y_new) > convergeThresh) then
               if (abs(y_prev - y_new) .gt. maxabsdiff) then
                  maxabsdiff = abs(y_prev - y_new)
                  imaxabsdiff = j
               end if
               done = .false.
            end if

            ! Updates mean absolute differences against observed values
            if (verbose) then
               icount = icount + 1
               y_est = pnew_est(j)
               y_ana = pnew_ana(j)
               y_obs = this%obs(j)
               diff = y_est - y_obs
               mad_est = mad_est + abs(diff)
               diff = y_ana - y_obs
               mad_ana = mad_ana + abs(diff)

               ! A crude check for unusually high normalized deviations.
               normdev = (y_obs - y_ana)*(y_obs - y_ana)/this%sigmaOSqr(j)
               if (normdev .gt. 9) then
                  num_high_dev = num_high_dev + 1
               end if

            end if

            pprev_est(j) = pnew_est(j)
            pprev_ana(j) = pnew_ana(j)

         end do ! j

         if (verbose) then
            if (icount .gt. 0) then
               mad_est = mad_est / real(icount)
               mad_ana = mad_ana / real(icount)
            end if
         end if

         ! !EMK TEST
         ! do ii = 1, nobs
         !    write(LIS_logunit,*) &
         !         '[INFO] ii,net,platform,obs, back, ana, est, dataDensity: ', &
         !         ii, &
         !         trim(this%net(ii)), ' ',&
         !         trim(this%platform(ii)), &
         !         ' ',this%obs(ii),&
         !         ' ',this%back(ii),&
         !         ' ',pnew_ana(ii),&
         !         ' ',pnew_est(ii),&
         !         ' ',1./invDataDensities(ii)
         ! end do

         if (done) exit ! No more iterations!

         if (verbose) then
            write(LIS_logunit,*) &
                 '[INFO] Bratseth scheme not converged yet after ',npasses, &
                 'iterations'
            write(LIS_logunit,*) &
                 '[INFO] Mean absolute difference against obs: ana: ', &
                 mad_ana,' est: ',mad_est
            write(LIS_logunit,*) &
                 '[INFO] Max abs change ',maxabsdiff,' at i = ', imaxabsdiff
            write(LIS_logunit,*) &
                 '[INFO] net,platform,obs, back, ana, est, dataDensity: ', &
                 trim(this%net(imaxabsdiff)), ' ',&
                 trim(this%platform(imaxabsdiff)), &
                 ' ',this%obs(imaxabsdiff),&
                 ' ',this%back(imaxabsdiff),&
                 ' ',pnew_ana(imaxabsdiff),&
                 ' ',pnew_est(imaxabsdiff),&
                 ' ',1./invDataDensities(imaxabsdiff)
           write(LIS_logunit,*) &
                '[INFO] Normalized deviation: ',normdev
           write(LIS_logunit,*) &
                '[INFO] Number of large normalized deviations: ',num_high_dev
           write(LIS_logunit,*) &
                '-------------------------------------------------------------'
        end if ! verbose

#if (defined SPMD)
         call MPI_Barrier(LIS_MPI_COMM, ierr)
         call handle_mpi_error(ierr, &
              'MPI_Barrier call in calc_obsAnalysis')
         t2 = MPI_Wtime()
         if (verbose) then
            write(LIS_logunit,*) &
                 '[INFO] Elapsed time for this iteration is ',t2 - t1, &
                 ' seconds'
         end if
#endif

	! EMK...Escape if too many iterations.  Since this is intended for
        ! production, we will allow the program to continue instead of
        ! aborting.
         !if (npasses .eq. 100) then
         if (npasses .eq. 5) then ! For Ops

            write(LIS_logunit,*) &
                 '[WARN] Bratseth failed to converge after ',npasses, &
                 ' iterations!'
            write(LIS_logunit,*) &
                 '[WARN] Will stop iterating'
            flush(LIS_logunit)
            exit
         end if
      end do ! Iterate until convergence

      if (verbose) then
         if (done) then
            write(LIS_logunit,*) &
              '[INFO] Bratseth analysis converged after ',npasses, &
              ' iterations'
         endif
         write(LIS_logunit,*) &
              '[INFO] Mean absolute difference against obs: ana: ',mad_ana, &
              ' est: ',mad_est
      end if

#if (defined SPMD)
      call MPI_Barrier(LIS_MPI_COMM, ierr)
      call handle_mpi_error(ierr, &
           'MPI_Barrier call in calc_obsAnalysis')
      t2 = MPI_Wtime()
      if (verbose) then
         write(LIS_logunit,*) &
              '[INFO] Elapsed time for final iteration is ',t2 - t1,' seconds'
         write(LIS_logunit,*) &
              '[INFO] Total elapsed time for all iterations is ',t2-t0, &
              'seconds'
         write(LIS_logunit,*) &
              '---------------------------------------------------------------'
      end if
#endif

      if (verbose) then
         icount = 0
         do j = 1, nobs
            if (this%qc(j) .eq. QC_REJECT) cycle

            if (present(skip)) then
               if (skip(j)) cycle
            end if

            icount = icount + 1
         end do
         write(LIS_logunit,*) '[INFO] ',icount,' obs used in this analysis'
      end if

      ! Collect the observed, background, and analysis values.
      if (agrmet_struc(nest)%oba_switch .eq. 1) then
         good_obs = USAF_countGoodObs(this)
         varOBA = createOBA(nest, maxobs=good_obs)
         do j = 1, nobs
            if (this%qc(j) .eq. QC_REJECT) cycle
            if (present(skip)) then
               if (skip(j)) cycle
            end if
            O = this%obs(j)
            B = this%back(j) ! Reusing B variable here
            A = pprev_ana(j)
            call assignOBA(varOBA, &
                 this%net(j), this%platform(j), &
                 this%lat(j), this%lon(j), &
                 O, B, A)
         end do ! j
      end if

      ! Clean up
      call destroy_hash2d(LIS_rc%gnc(nest),LIS_rc%gnr(nest),hash2d)
      deallocate(hash2d)
      deallocate(sumObsEstimates_pet)
      deallocate(pprev_est)
      deallocate(pprev_ana)
      deallocate(pnew_est)
      deallocate(pnew_est_pet)
      deallocate(pnew_ana)
      deallocate(pnew_ana_pet)

   end subroutine calc_obsAnalysis

   !---------------------------------------------------------------------------
   ! Perform Bratseth analysis at grid points.  Assumes (1) the mrgp array
   ! contains the transformed background first guess; (2) the Bratseth scheme
   ! was already run at the observation points; and (3) the summed observation
   ! estimates and number of passes from that operation is provided (as
   ! sumObsEstimates and npasses, respectively).
   !
   ! The interpolation from observation points to grid points is done in
   ! a single pass, similar to Daley (1991) or Kalnay (2003).  This greatly
   ! saves time compared to the original Bratseth (1986) or Sashegyi et al
   ! (1993) approaches, where ob-to-grid interpolation was done as each
   ! analysis pass was performed at the observation points.
   !
   ! This implementation is parallelized by farming work out to the LIS ESMF
   ! PETs.  To improve performance, a 2D hash table is constructed to group
   ! observations by LIS grid box; this allows the code to avoid ob-to-ob
   ! comparisons that are obviously too far away to be correlated.
   !
   ! NOTE:  Bratseth values are not interpolated to water points.

   subroutine calc_gridAnalysis(this,nest,sigmaBSqr,nobs,invDataDensities,&
        sumObsEstimates,npasses,back,max_dist,backErrScaleLength,mrgp)

      ! Imports
      use LIS_coreMod, only: LIS_rc, LIS_domain, &
           LIS_ews_halo_ind, LIS_ewe_halo_ind, &
           LIS_nss_halo_ind, LIS_nse_halo_ind, LIS_localPet
      use LIS_logMod, only : LIS_logunit, LIS_endrun
      use LIS_LMLCMod, only: LIS_LMLC
      use LIS_mpiMod

      ! Defaults
      implicit none

      ! Arguments
      type(USAF_ObsData),intent(inout) :: this
      integer, intent(in) :: nest
      real, intent(in) :: sigmaBSqr
      integer, intent(in) :: nobs
      real, intent(in) :: invDataDensities(nobs)
      real, intent(in) :: sumObsEstimates(nobs)
      integer, intent(in) :: npasses
      real, intent(in) :: back(LIS_rc%gnc(nest), LIS_rc%gnr(nest))
      real, intent(in) :: max_dist
      real, intent(in) :: backErrScaleLength
      real, intent(inout) :: mrgp(LIS_rc%lnc(nest),LIS_rc%lnr(nest))

      ! Local variables
      type(hash_list), allocatable, target :: hash2d(:,:)
      integer :: cmax,rmax
      double precision :: t1, t2
      real ::  locallat, locallon, tmp_mrgp
      real :: dist, weight
      integer :: pet, pet_incr
      integer :: r, c, j, job, ierr
      integer :: lr,ur,lc1,rc1,lc2,rc2
      integer :: nobs_neighbors
      integer, allocatable :: jobs_neighbors_vector(:)
      real, allocatable :: mrgp_1d_pet(:), mrgp_1d(:)
      real :: back1(1)
      integer :: r_local, c_local
      integer :: gindex

      ! Sanity checks
      if (nobs .eq. 0) return

      if (nobs .ne. this%nobs) then
         write(LIS_logunit,*)'[ERR] Array size mismatch in calc_gridAnalysis!'
         write(LIS_logunit,*)'[ERR] nobs, this%nobs = ',nobs, this%nobs
         call LIS_endrun()
      end if

#if (defined SPMD)
      call MPI_Barrier(LIS_MPI_COMM, ierr)
      call handle_mpi_error(ierr, &
              'MPI_Barrier call in calc_gridAnalysis')
      t1 = MPI_Wtime()
#endif

      ! Here we create a 2d hash table storing the index values of each ob
      ! in linked lists for each LIS grid box.  This can help us screen
      ! out obviously unnecessary ob comparisons later.
      call build_hash2d(this,nest,cmax,rmax,hash2d)

      write(LIS_logunit,*) &
           '[INFO] Calculating analysis at grid points...'

      allocate(mrgp_1d_pet(LIS_rc%gnc(nest)*LIS_rc%gnr(nest)))
      mrgp_1d_pet(:) = 0

      pet = -1
      pet_incr = 1

      ! Now calculate the analysis at each grid point.
      do r = 1,LIS_rc%gnr(nest)
         do c = 1,LIS_rc%gnc(nest)

            ! Skip if water point
            ! FIXME -- Allow use over water?
            if (LIS_LMLC(nest)%glandmask(c,r) .le. 0) cycle

            ! See which PET is responsible for this grid box.
            call update_pet(pet,pet_incr)
            if (pet .ne. LIS_localPet) cycle

            gindex = c+(r-1)*LIS_rc%gnc(nest)

            locallat = LIS_domain(nest)%glat(gindex)
            locallon = LIS_domain(nest)%glon(gindex)

            ! Find neighbors with positive correlation in background field.
            call find_gridpt_neighbors(c,r,nest,max_dist, &
                 lr,ur,lc1,rc1,lc2,rc2)

            ! Get list of obs from neighboring grid boxes
            call get_neighbor_obs(LIS_rc%gnc(nest),LIS_rc%gnr(nest),hash2d,&
                 lr,ur,lc1,rc1,lc2,rc2,&
                 nest, nobs_neighbors, jobs_neighbors_vector)

            back1(1) = back(c,r)
            tmp_mrgp = back1(1)

            if (nobs_neighbors .eq. 0) then
               mrgp_1d_pet(gindex) = tmp_mrgp
               cycle
            end if

            do j = 1, nobs_neighbors

               job = jobs_neighbors_vector(j)

               ! Skip bad observations.
               if (this%qc(job) .eq. QC_REJECT) cycle

               dist = great_circle_distance(locallat,locallon, &
                    this%lat(job), this%lon(job))
               if (dist .gt. max_dist) cycle

               weight = &
                    backErrCov(sigmaBSqr,dist,backErrScaleLength) &
                    * invDataDensities(job)

               tmp_mrgp = tmp_mrgp + &
                    (weight * ((npasses*this%obs(job)) - sumObsEstimates(job)))
            end do ! j

            mrgp_1d_pet(gindex) = tmp_mrgp

            deallocate(jobs_neighbors_vector)

         end do ! c
      end do ! r

      ! Clean up
      call destroy_hash2d(LIS_rc%gnc(nest),LIS_rc%gnr(nest),hash2d)
      deallocate(hash2d)

#if (defined SPMD)
      allocate(mrgp_1d(LIS_rc%gnc(nest)*LIS_rc%gnr(nest)))
      mrgp_1d(:) = 0
      call MPI_Barrier(LIS_MPI_COMM, ierr)
      call handle_mpi_error(ierr, &
           'MPI_Barrier call in calc_gridAnalysis')
      call MPI_ALLREDUCE(mrgp_1d_pet,mrgp_1d, &
           LIS_rc%gnc(nest)*LIS_rc%gnr(nest),MPI_REAL, &
           MPI_SUM, LIS_MPI_COMM, ierr)
      call handle_mpi_error(ierr, &
           'MPI_ALLREDUCE call in calc_gridAnalysis')
      deallocate(mrgp_1d_pet)
#endif

      ! Now copy to local 2D array
      do r = LIS_nss_halo_ind(nest,LIS_localPet+1), &
           LIS_nse_halo_ind(nest,LIS_localPet+1)
         r_local = r - LIS_nss_halo_ind(nest,LIS_localPet+1) + 1

         do c = LIS_ews_halo_ind(nest,LIS_localPet+1), &
              LIS_ewe_halo_ind(nest,LIS_localPet+1)

            ! EMK...Make sure this is a land point
            if (LIS_LMLC(nest)%glandmask(c,r) .le. 0) cycle

            c_local = c - LIS_ews_halo_ind(nest,LIS_localPet+1) + 1

            gindex = c+(r-1)*LIS_rc%gnc(nest)

            mrgp(c_local,r_local) = mrgp_1d(gindex)

         end do ! c
      end do ! r
      deallocate(mrgp_1d)

#if (defined SPMD)
      call MPI_Barrier(LIS_MPI_COMM, ierr)
      call handle_mpi_error(ierr, &
           'MPI_Barrier call in calc_gridAnalysis')
      t2 = MPI_Wtime()
      write(LIS_logunit,*)'[INFO] Elapsed time for grid analysis is ', &
           t2 - t1,' seconds'
#endif

   end subroutine calc_gridAnalysis

   !---------------------------------------------------------------------------
   ! Observation error covariance function.
   real function obsErrCov(sigmaOSqr,oErrScaleLength,dist)
      implicit none
      real, intent(in) :: sigmaOSqr
      real, intent(in) :: oErrScaleLength ! in meters
      real, intent(in) :: dist ! in meters
      obsErrCov = sigmaOSqr*obsErrCorr(oErrScaleLength,dist)
   end function obsErrCov

   !---------------------------------------------------------------------------
   ! Observation error correlation function.  Currently Gaussian.
   real function obsErrCorr(oErrScaleLength,dist)
      implicit none
      real, intent(in) :: oErrScaleLength ! in meters
      real, intent(in) :: dist ! in meters
      real :: invOErrScaleLength
      invOErrScaleLength = 1. / oErrScaleLength
      obsErrCorr = exp(-1*dist*dist*invOErrScaleLength*InvOErrScaleLength)
   end function obsErrCorr

   !---------------------------------------------------------------------------
   ! Background error covariance function.
   real function backErrCov(sigmaBSqr,dist,scale_length)
      implicit none
      real, intent(in) :: sigmaBSqr
      real, intent(in) :: dist ! in meters
      real, intent(in) :: scale_length
      backErrCov = sigmaBSqr*backErrCorr(dist,scale_length)
   end function backErrCov

   !---------------------------------------------------------------------------
   ! Background error correlation function.  Currently Gaussian.
   real function backErrCorr(dist,scale_length)
      implicit none
      real, intent(in) :: dist ! in meters
      real, intent(in) :: scale_length
      real :: inv_scale_length
      inv_scale_length = 1./scale_length
      backErrCorr = exp(-1*dist*dist*inv_scale_length*inv_scale_length)
   end function backErrCorr

   !---------------------------------------------------------------------------
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

   !---------------------------------------------------------------------------
   ! Read precipitation fields from GFS or GALWEM, split into 3-hr
   ! accumulation, and interpolate to LIS global domain array.
   subroutine fldbld_precip_nwp(nest,findex,julhr,src,fc_hr, &
        nwp_precip_glb,rc)

      ! Imports
      use AGRMET_forcingMod, only : agrmet_struc
#if (defined USE_GRIBAPI)
      use grib_api
#else
      use LIS_coreMod,       only : LIS_masterproc
#endif
      use LIS_coreMod,       only : LIS_rc
      use LIS_logMod,        only : LIS_logunit, LIS_abort, LIS_alert, &
           LIS_verify, LIS_endrun
      use LIS_timeMgrMod,    only : LIS_julhr_date
      use LIS_historyMod,    only : LIS_gather_2d_local_to_global

      ! Defaults
      implicit none

      ! Arguments
      integer, intent(in) :: nest
      integer, intent(in) :: findex
      integer, intent(in) :: julhr
      character(len=6), intent(in) :: src ! "GFS" or "GALWEM"
      integer, intent(inout)   :: fc_hr
      real, intent(out)   :: nwp_precip_glb(LIS_rc%gnc(nest), LIS_rc%gnr(nest))
      integer, intent(out) :: rc

      ! Local variables
      logical :: found, found2
      integer :: yr1, mo1, da1, hr1
      integer :: file_julhr
      integer :: getsixhr
      integer :: yr_2d
      !character(len=120) :: gribfile, gribfile2
      character(len=255) :: gribfile, gribfile2
      integer :: center
      real :: gridres
      integer :: Ni, Nj, ifguess, jfguess
      real, allocatable :: fg_prec1(:,:), fg_prec2(:,:), fg_prec(:,:)
      real, allocatable :: fg_data(:,:)
      integer :: rc2
      logical :: first_time

      external :: getAVNfilename
      external :: AGRMET_getGALWEMfilename
      external :: AGRMET_fg2lis_precip
      external :: interp_galwem_first_guess

      rc = 0

      found = .false.
      found2 = .false.
      call LIS_julhr_date(julhr,yr1,mo1,da1,hr1)
      file_julhr = julhr
      first_time = .true.

      ! Need to process current, and sometimes previous, accumuluations.
      ! Search for an analysis of forecast file for up to 24 hours with the
      ! needed valid time.
      do while( ((.not.found) .or. (.not.found2)) .and. (fc_hr .le. 24))

         ! Make sure start with previous cycle of GRIB data when using
         ! GALWEM.  But for GFS, we can check the current cycle for
         ! legacy reasons.
         if ( (.not. first_time) .or. &
              (first_time .and. trim(src) == "GALWEM" .and. fc_hr .le. 6)) then
            fc_hr = fc_hr + 6
            if (fc_hr > 24) exit ! Give up
            file_julhr = file_julhr - 6
            call LIS_julhr_date(file_julhr,yr1,mo1,da1,hr1)
         end if
         first_time = .false.
         found = .false.

         ! Some GFS files have 6-hr accumulations, which require also
         ! reading a 3-hr accumulation from a previous file and differencing
         ! (e.g., read 0-6 hr, read earlier 0-3 hr, and subtract to get 3-6
         ! hr).  See if we need the 6-hour forecast file.
         if (mod(fc_hr,6) .eq. 0 .and. trim(src) == "GFS") then
            getsixhr=1      ! We will read a 6-hr file
            found2=.false.  ! Will search for separate 3-hr file
         else
            getsixhr=0      ! We will not read a 6-hr file
            found2=.true.   ! No need to search for separate 3-hr file
         endif

         yr_2d = mod(yr1,100)
         if (yr_2d.eq.0) yr_2d = 100
         if (src .eq. "GFS") then
            ! EMK...Added support for new GFS filename convention
            call getAVNfilename(gribfile, agrmet_struc(nest)%agrmetdir,&
                 agrmet_struc(nest)%gfsdir, agrmet_struc(nest)%use_timestamp,&
                 agrmet_struc(nest)%gfs_timestamp, &
                 agrmet_struc(nest)%gfs_filename_version, &
                 yr1, mo1, da1, hr1, fc_hr)
            if (getsixhr.eq.1) then
               call getAVNfilename(gribfile2, agrmet_struc(nest)%agrmetdir,&
                    agrmet_struc(nest)%gfsdir,&
                    agrmet_struc(nest)%use_timestamp,&
                    agrmet_struc(nest)%gfs_timestamp,&
                    agrmet_struc(nest)%gfs_filename_version, &
                    yr1, mo1, da1, hr1, fc_hr-3)
            endif
         else if (src .eq. "GALWEM") then
            call AGRMET_getGALWEMfilename(gribfile, &
                 agrmet_struc(nest)%agrmetdir,&
                 agrmet_struc(nest)%galwemdir,&
                 agrmet_struc(nest)%use_timestamp,&
                 agrmet_struc(nest)%galwem_res, &
                 yr1,mo1,da1,hr1,fc_hr)
         end if

         ! Check for valid time in first GRIB file
         call check_grib_file(gribfile,yr1,mo1,da1,hr1,found, &
              center, Ni, Nj, gridres)
         if (.not. found) then
            write(LIS_logunit,*) &
                 '[WARN] Problem finding valid time in GRIB file: ', &
                 trim(gribfile)
            cycle ! Roll back
         end if

         ! Check for valid time in second GRIB file if needed
         if (getsixhr .eq. 1) then
            call check_grib_file(gribfile2,yr1,mo1,da1,hr1,found2, &
                 center, Ni, Nj, gridres)
            if (.not. found2) then
               write(LIS_logunit,*) &
                    '[WARN] Problem finding valid time in GRIB file: ', &
                    trim(gribfile2)
               cycle ! Roll back
            end if
         end if

         ! At this point, the file(s) have the correct time.  Try extracting
         ! the precipitation data.
         write(LIS_logunit,*)'[INFO] FIRST GUESS DATA IS ON A ', gridres,&
              ' DEGREE LAT/LON GRID'
         ifguess = Ni
         jfguess = Nj
         if (center .eq. 7) then
            write(LIS_logunit,*)'[INFO] FIRST GUESS DATA IS FROM GFS MODEL'
         elseif (center .eq. 57) then
            write(LIS_logunit,*) &
                 '[INFO] FIRST GUESS DATA IS FROM UK UM (GALWEM) MODEL'
         elseif (center .eq. 58) then
            write(LIS_logunit,*)'[INFO] FIRST GUESS DATA IS FROM NOGAPS MODEL'
         end if

         if (src .eq. 'GFS') then
            call fldbld_read_precip_gfs( gribfile, ifguess, jfguess,&
                 fg_prec1, rc2 )
            if (rc2 .ne. 0) then
               write(LIS_logunit,*) &
                    '[WARN] Problem reading precipitation from GRIB file: ', &
                    trim(gribfile)
               found = .false.
               found2 = .false.
               if (allocated(fg_prec1)) deallocate(fg_prec1)
               cycle ! Roll back
            end if
            if (getsixhr.eq.1) then
               call fldbld_read_precip_gfs( gribfile2, ifguess, &
                    jfguess, fg_prec2, rc2 )
               if (rc2 .ne. 0) then
                  write(LIS_logunit,*) &
                       '[WARN] Problem reading precipitation from GRIB file: '&
                       , trim(gribfile2)
                  found = .false.
                  found2 = .false.
                  if (allocated(fg_prec1)) deallocate(fg_prec1)
                  if (allocated(fg_prec2)) deallocate(fg_prec2)
                  cycle ! Roll back
               end if
            end if
         else if (src .eq. 'GALWEM') then
            call fldbld_read_precip_galwem(gribfile, ifguess, jfguess, fc_hr, &
                 fg_prec1, rc2)
            if (rc2 .ne. 0) then
               write(LIS_logunit,*) &
                    '[WARN] Problem reading precipitation from GRIB file: ', &
                    trim(gribfile)
               found = .false.
               found2 = .false.
               if (allocated(fg_prec1)) deallocate(fg_prec1)
               cycle ! Roll back
            end if
         end if

      end do ! do while

      ! Handle missing GRIB data
      if ( (.not. found) .or. (.not. found2) ) then
         if (src .eq. "GFS") then
            write(LIS_logunit,*) &
                 '[ERR] ** GFS Precipitation data not available **'
         else if (src .eq. "GALWEM") then
            write(LIS_logunit,*) &
                 '[ERR] ** GALWEM Precipitation data not available **'
         end if
         rc = 1
         return
      end if

      ! Log which grib file(s) was/were ultimately selected
      write(LIS_logunit,*) &
           '[INFO] Using NWP precipitation from ',trim(gribfile)
      if (getsixhr.eq.1) then
         write(LIS_logunit,*) &
              '[INFO] Also using NWP precipitation from ',trim(gribfile2)
         write(LIS_logunit,*) &
              '[INFO] Will difference two files to get 3-hr accumulation'
      end if

      allocate ( fg_prec (ifguess, jfguess) )
      if (getsixhr.eq.1) then
         fg_prec = fg_prec1 - fg_prec2
      else
         fg_prec = fg_prec1
      endif

      ! Sometimes subtraction of 3-hr precip from 6-hr causes slightly
      ! negative values.  Correct this.
      where (fg_prec .lt. 0)
         fg_prec=0
      endwhere

      ! Interpolate to the LIS grid
      allocate ( fg_data (LIS_rc%lnc(nest), LIS_rc%lnr(nest)) )
      if ( src .eq. 'GFS' ) then
         call AGRMET_fg2lis_precip(nest, findex, ifguess, jfguess, &
                                   fg_prec, fg_data)
      else if ( src .eq. 'GALWEM' ) then
         call interp_galwem_first_guess(nest, ifguess, jfguess, .true., &
                                        fg_prec, fg_data)
      endif
      call LIS_gather_2d_local_to_global(nest, fg_data, nwp_precip_glb)

      ! Clean up
      deallocate ( fg_prec )
      deallocate ( fg_prec1 )
      deallocate ( fg_data )
      if (getsixhr.eq.1) deallocate ( fg_prec2 )
      rc = 0

   end subroutine fldbld_precip_nwp

   !---------------------------------------------------------------------------
   ! Checks GFS or GALWEM GRIB file and retrieves originating center, date,
   ! time, grid resolution, and grid dimensions.
   subroutine check_grib_file(gribfile,yr1,mo1,da1,hr1,found, &
        center, Ni, Nj, gridres)

      ! Imports
#if (defined USE_GRIBAPI)
      use grib_api
#endif
      use LIS_coreMod, only: LIS_masterproc
      use LIS_logMod,  only : LIS_logunit, LIS_abort, LIS_alert, &
           LIS_verify, LIS_endrun
      use LIS_mpiMod

      ! Defaults
      implicit none

      ! Arguments
      !character(len=120), intent(in) :: gribfile
      character(len=255), intent(in) :: gribfile
      integer, intent(in) :: yr1, mo1, da1, hr1
      logical, intent(out) :: found
      integer, intent(out) :: center
      integer, intent(out) :: Ni, Nj
      real, intent(out) :: gridres

      ! Local variables
      integer :: ftn, ierr, igrib
      integer :: dataDate, dataTime
      character(len=100) :: gtype
      logical :: found_inq
#if (!defined USE_GRIBAPI)
      character(len=255) :: message(20)
#endif
      found = .false.
      ! Dummy values, replaced by GRIB data contents below
      center = 0
      Ni = 0
      Nj = 0
      gridres = 0

      ! EMK...Before using ECCODES/GRIB_API, see if the GRIB file exists
      ! using a simple inquire statement.  This avoids ECCODES/GRIB_API
      ! writing error messages to stdout/stderr, which may lead to runtime
      ! problems.
      inquire(file=trim(gribfile),exist=found_inq)
      if (.not. found_inq) then
         write(LIS_logunit,*)'[WARN] Cannot find file '//trim(gribfile)
         return
      end if

#if (defined USE_GRIBAPI)
      call grib_open_file(ftn,trim(gribfile),'r',ierr)
      if ( ierr .ne. 0 ) then
         write(LIS_logunit,*) '[WARN] Failed to open - ', trim(gribfile)
         return
      end if

      call grib_new_from_file(ftn,igrib,ierr)
      if ( ierr .ne. 0 ) then
         write(LIS_logunit,*) '[WARN] failed to read - '//trim(gribfile)
         call grib_release(igrib,ierr)
         call grib_close_file(ftn)
         return
      endif

      call grib_get(igrib,'centre',center,ierr)
      if ( ierr .ne. 0 ) then
         write(LIS_logunit,*) '[WARN] error in grib_get: centre in ' // &
              'check_grib_file'
         call grib_release(igrib,ierr)
         call grib_close_file(ftn)
         return
      endif

      call grib_get(igrib,'gridType',gtype,ierr)
      if ( ierr .ne. 0 ) then
         write(LIS_logunit,*) '[WARN] error in grid_get: gridtype in ' // &
              'check_grib_file'
         call grib_release(igrib,ierr)
         call grib_close_file(ftn)
         return
      endif
      if ( gtype .ne. "regular_ll" ) then
         write(LIS_logunit,*) &
              '[WARN] GRIB data not on regular lat-lon grid!'
         call grib_release(igrib,ierr)
         call grib_close_file(ftn)
         return
      endif

      call grib_get(igrib,'Ni',Ni,ierr)
      if ( ierr .ne. 0 ) then
         write(LIS_logunit,*) '[WARN] error in grid_get:Ni in ' // &
              'check_grib_file'
         call grib_release(igrib,ierr)
         call grib_close_file(ftn)
         return
      endif

      call grib_get(igrib,'Nj',Nj,ierr)
      if ( ierr .ne. 0 ) then
         write(LIS_logunit,*) '[WARN] error in grid_get:Nj in ' // &
              'check_grib_file'
         call grib_release(igrib,ierr)
         call grib_close_file(ftn)
         return
      endif

      call grib_get(igrib,'jDirectionIncrementInDegrees',gridres,ierr)
      if ( ierr .ne. 0 ) then
         write(LIS_logunit,*) &
              '[WARN] error in grid_get:jDirectionIncrementInDegrees in ' // &
              'check_grib_file'
         call grib_release(igrib,ierr)
         call grib_close_file(ftn)
         return
      endif

      call grib_get(igrib,'dataDate',dataDate,ierr)
      if ( ierr .ne. 0 ) then
         write(LIS_logunit,*) '[WARN] error in grid_get:dataDate in ' // &
              'check_grib_file'
         call grib_release(igrib,ierr)
         call grib_close_file(ftn)
         return
      endif

      call grib_get(igrib,'dataTime',dataTime,ierr)
      if ( ierr .ne. 0 ) then
         write(LIS_logunit,*) '[WARN] error in grid_get:dataTime in ' // &
              'check_grib_file'
         call grib_release(igrib,ierr)
         call grib_close_file(ftn)
         return
      endif

      if ( yr1*10000+mo1*100+da1 .ne. dataDate) then
         write(LIS_logunit,*) '[WARN] Cannot find valid date in GRIB file!'
         call grib_release(igrib,ierr)
         call grib_close_file(ftn)
         return
      end if

      if (hr1*100 .ne. dataTime ) then
         write(LIS_logunit,*) '[WARN] Cannot find valid time in GRIB file!'
         call grib_release(igrib,ierr)
         call grib_close_file(ftn)
         return
      endif

      call grib_release(igrib,ierr)
      if ( ierr .ne. 0 ) then
         write(LIS_logunit,*) &
              '[WARN] error in grid_release in ' // &
              'check_grib_file'
         return
      end if

      ! At this point, we declare victory
      call grib_close_file(ftn)
      found = .TRUE.

#else
      write(LIS_logunit,*) &
           '[ERR] check_grib_file requires GRIB-API or ECCODES library'
      write(LIS_logunit,*) '[ERR] please recompile LIS'
      flush(LIS_logunit)
      message(:) = ''
      message(1) = '[ERR] Program:  LIS'
      message(2) = '  Routine: check_grib_file.'
      message(3) = '  LIS was not compiled with GRIBAPI or ECCODES support!'
#if (defined SPMD)
      call MPI_Barrier(LIS_MPI_COMM, ierr)
      call handle_mpi_error(ierr, &
           'MPI_Barrier call in check_grib_file')
#endif
      if(LIS_masterproc) then
         call LIS_alert( 'LIS.check_grib_file', 1, &
              message )
         call LIS_abort( message)
      else
         call sleep(10) ! Make sure LIS_masterproc finishes LIS_abort
         call LIS_endrun()
      endif

      call LIS_endrun()
#endif

   end subroutine check_grib_file

   !---------------------------------------------------------------------------
   ! Creates "superobs" out of close observations.  Each close observation is
   ! first checked for unacceptable deviation from the mean of the close
   ! obs, and rejected if deviation is too high.  Superobs are considered
   ! "close" if they are in the same LIS grid box.  Based on Lespinas et al
   ! (2015).
   !
   ! This implementation is parallelized by distributing work among the
   ! different LIS ESMF PETs.
   subroutine USAF_superstatQC(this,nest,new_name,network,silent_rejects)

      ! Imports
      use LIS_coreMod, only: LIS_domain, LIS_rc, LIS_localPet
      use LIS_logMod, only: LIS_logunit, LIS_endrun
      use LIS_mpiMod

      ! Defaults
      implicit none

      ! Arguments
      type(USAF_ObsData), intent(inout) :: this
      integer,intent(in) :: nest
      character(len=32), intent(in) :: new_name
      character(len=*), optional :: network
      logical,optional,intent(in) :: silent_rejects

      ! Local variables
      integer :: nobs
      integer :: num_rejected_obs, num_merged_obs, num_superobs
      real :: dlat, dlon
      integer :: c,r,j
      real :: ctrlat, ctrlon
      character(len=32) :: net_new
      character(len=32) :: platform_new
      integer, allocatable :: actions(:), actions_pet(:)
      real, allocatable :: superobs_pet(:),superlat_pet(:),superlon_pet(:)
      real, allocatable :: means(:)
      real, allocatable :: superobs(:),superlat(:),superlon(:)
      real, allocatable :: superSigmaOSqr(:), superSigmaOSqr_pet(:)
      real, allocatable :: superOErrScaleLength(:), superOErrScaleLength_pet(:)
      integer, allocatable :: superob_count(:), superob_count_pet(:)
      integer :: pet, pet_incr
      integer :: glbcr
      integer :: ierr, ipass, icount
      integer :: gindex
      real :: threshold
      logical :: found
      double precision :: t1, t2
      logical :: silent_rejects_local

      ! Sanity check
      nobs = this%nobs
      if (nobs .eq. 0) then
         write(LIS_logunit,*)&
              '[INFO] superstatQC found no observations to test'
         return
      end if

      silent_rejects_local = .false.
      if (present(silent_rejects)) then
         silent_rejects_local = silent_rejects
      end if

      ! Determine which type of observations will be "superobbed"
      net_new = trim(new_name)
      platform_new = trim(new_name)

!       if (present(network)) then
!          write(LIS_logunit,*)&
!               '[INFO] superstatQC working with ',nobs,' obs from ', &
!               trim(network),'...'
!       else
!          write(LIS_logunit,*)&
!               '[INFO] superstatQC working with ',nobs,' obs...'
!       end if

#if (defined SPMD)
      call MPI_Barrier(LIS_MPI_COMM, ierr)
      call handle_mpi_error(ierr, &
           'MPI_Barrier call in USAF_superstatQC')
      t1 = MPI_Wtime()
#endif

      allocate(actions_pet(nobs))
      actions_pet(:) = 0
      allocate(actions(nobs))
      actions(:) = 0

      glbcr = LIS_rc%gnc(nest)*LIS_rc%gnr(nest)
      allocate(superobs_pet(glbcr))
      superobs_pet(:) = 0
      allocate(superlat_pet(glbcr))
      superlat_pet(:) = 0
      allocate(superlon_pet(glbcr))
      superlon_pet(:) = 0
      allocate(superSigmaOSqr_pet(glbcr))
      superSigmaOSqr_pet(:) = 0
      allocate(superOErrScaleLength_pet(glbcr))
      superOErrScaleLength_pet(:) = 0
      allocate(superob_count_pet(glbcr))
      superob_count_pet(:) = 0

      allocate(means(glbcr))
      means(:) = 0
      allocate(superobs(glbcr))
      superobs(:) = 0
      allocate(superlat(glbcr))
      superlat(:) = 0
      allocate(superlon(glbcr))
      superlon(:) = 0
      allocate(superSigmaOSqr(glbcr))
      superSigmaOSqr(:) = 0
      allocate(superOErrScaleLength(glbcr))
      superOErrScaleLength(:) = 0
      allocate(superob_count(glbcr))
      superob_count(:) = 0

      dlat = LIS_domain(nest)%lisproj%dlat
      dlon = LIS_domain(nest)%lisproj%dlon

      ! First pass: Find all acceptable obs in each LIS grid box, and
      ! calculate average value.
      !
      ! Second pass:  Reject observations that deviate too much from
      ! local average.  Calculate superob value, location, and error variance.
      !
      ! Third pass:  Make sure obs that could not be merged (due to
      ! lack of valid neighbors) are marked for preservation.

      do ipass = 1,3

         pet = -1
         pet_incr = 1

         do j = 1, nobs

            ! See which MPI process is responsible for this ob
            call update_pet(pet,pet_incr)
            if (pet .ne. LIS_localPet) cycle

            ! Skip bad data
            if ( this%qc(j) .eq. QC_REJECT) cycle

            ! Screen by type
            if (present(network)) then
               if (trim(network) .eq. "SSMI") then
                  if (.not. is_ssmi(this%net(j))) cycle
               else if (trim(network) .eq. "GEOPRECIP") then
                  if (.not. is_geoprecip(this%net(j))) cycle
               else if (trim(network) .eq. "CMORPH") then
                  if (.not. is_cmorph(this%net(j))) cycle
               else if (trim(network) .eq. "IMERG") then
                  if (.not. is_imerg(this%net(j))) cycle
               end if
            else ! Gauges
               if (.not. is_gauge(this%net(j))) cycle
            end if

            ! Now see which LIS grid box this is in.  First, handle latitude.
            found = .false.
            do r = 1, LIS_rc%gnr(nest)
               gindex = 1 + (r-1)*LIS_rc%gnc(nest)
               ctrlat = LIS_domain(nest)%glat(gindex)
               if (r .eq. 1) then
                  if (this%lat(j) .lt. (ctrlat - (0.5*dlat))) cycle
               end if
               if (this%lat(j) .ge. (ctrlat + (0.5*dlat))) cycle
               found = .true.
               exit
            end do ! r
            if (.not. found) cycle

            ! Now, longitude
            found = .false.
            do c = 1, LIS_rc%gnc(nest)
               gindex = c + (r-1)*LIS_rc%gnc(nest)
               ctrlon = LIS_domain(nest)%glon(gindex)
               if (c .eq. 1) then
                  if (this%lon(j) .lt. (ctrlon - (0.5*dlon))) cycle
               end if
               if (this%lon(j) .ge. (ctrlon + (0.5*dlon))) cycle
               found = .true.
               exit
            end do ! c
            if (.not. found) cycle

            ! Add contribution to local mean.
            if (ipass .eq. 1) then
               superobs_pet(gindex) = &
                    superobs_pet(gindex) + this%obs(j)
               superob_count_pet(gindex) = &
                    superob_count_pet(gindex) + 1
            end if ! ipass .eq. 1

            ! Reject obs that deviate too much from local mean, and collect
            ! contributions from the remainder for a superob.
            if (ipass .eq. 2) then
               if (means(gindex) .eq. MISSING) cycle

               icount = superob_count(gindex)
               threshold = 3 * this%sigmaOSqr(j) * &
                    sqrt(real(icount) / real(icount-1))

               if (abs(means(gindex) - this%obs(j)) .gt. threshold) then
                  actions_pet(j) = -1 ! Reject
               else
                  actions_pet(j) =  1 ! Consider for superob
                  superobs_pet(gindex) = &
                       superobs_pet(gindex) + this%obs(j)
                  superlat_pet(gindex) = &
                       superlat_pet(gindex) + this%lat(j)
                  superlon_pet(gindex) = &
                       superlon_pet(gindex) + this%lon(j)
                  superSigmaOSqr_pet(gindex) = &
                       superSigmaOSqr_pet(gindex) + this%sigmaOSqr(j)
                  superOErrScaleLength_pet(gindex) = &
                       superOErrScaleLength_pet(gindex) + &
                       this%oErrScaleLength(j)
                  superob_count_pet(gindex) = &
                       superob_count_pet(gindex) + 1
               end if ! Threshold check
            end if ! ipass .eq. 2

            ! If only one good observation is in this grid box, make sure
            ! it is preserved as is.
            if (ipass .eq. 3) then
               if (superobs(gindex) .eq. MISSING) then
                  if (actions_pet(j) .eq. 1) then
                     actions_pet(j) = 0
                  end if
               end if
            end if ! ipass .eq. 3

         end do ! j

         ! Collect contributions to local means
         if (ipass .eq. 1) then
            means(:) = 0
#if (defined SPMD)
            call MPI_Allreduce(superobs_pet,means,glbcr,MPI_REAL, MPI_SUM, &
                 LIS_MPI_COMM, ierr)
            call handle_mpi_error(ierr, &
                 'MPI_Allreduce call in USAF_superstatQC')
#else
            do j = 1, glbcr
               means(j) = sum(superobs_pet)
            end do
#endif
            superobs_pet(:) = 0
            superob_count(:) = 0
#if (defined SPMD)
            call MPI_Allreduce(superob_count_pet,superob_count,glbcr,&
                 MPI_INTEGER, MPI_SUM, LIS_MPI_COMM, ierr)
            call handle_mpi_error(ierr, &
                 'MPI_Allreduce call in USAF_superstatQC')
#else
            do j = 1, glbcr
               superob_count(j) = sum(superob_count_pet)
            end do
#endif
            superob_count_pet(:) = 0

            do r = 1, LIS_rc%gnr(nest)
               do c = 1, LIS_rc%gnc(nest)
                  gindex = c + (r-1)*LIS_rc%gnc(nest)
                  if (superob_count(gindex) .gt. 2) then
                     means(gindex) = means(gindex) / superob_count(gindex)
                  else ! Superob either not possible or not necessary
                     means(gindex) = MISSING
                  end if
               end do ! c
            end do ! r
         end if ! ipass.eq.1

         ! Collect contributions to local superob
         if (ipass .eq. 2) then

#if (defined SPMD)
            call MPI_Barrier(LIS_MPI_COMM, ierr)
            call handle_mpi_error(ierr, &
                 'MPI_Barrier call in USAF_superstatQC')
            superob_count(:) = 0
            call MPI_Allreduce(superob_count_pet, superob_count, glbcr,&
                 MPI_INTEGER, MPI_SUM, LIS_MPI_COMM, ierr)
            call handle_mpi_error(ierr, &
                 'MPI_Allreduce call in USAF_superstatQC')
            superob_count_pet(:) = 0

            superobs(:) = 0
            call MPI_Allreduce(superobs_pet, superobs, glbcr, MPI_REAL, &
                 MPI_SUM, LIS_MPI_COMM, ierr)
            call handle_mpi_error(ierr, &
                 'MPI_Allreduce call in USAF_superstatQC')
            superobs_pet(:) = 0

            superlat(:) = 0
            call MPI_Allreduce(superlat_pet, superlat, glbcr, MPI_REAL, &
                 MPI_SUM, LIS_MPI_COMM, ierr)
            call handle_mpi_error(ierr, &
                 'MPI_Allreduce call in USAF_superstatQC')
            superlat_pet(:) = 0

            superlon(:) = 0
            call MPI_Allreduce(superlon_pet, superlon, glbcr, MPI_REAL, &
                 MPI_SUM, LIS_MPI_COMM, ierr)
            call handle_mpi_error(ierr, &
                 'MPI_Allreduce call in USAF_superstatQC')
            superlon_pet(:) = 0

            superSigmaOSqr(:) = 0
            call MPI_Allreduce(superSigmaOSqr_pet, superSigmaOSqr, glbcr, &
                 MPI_REAL, MPI_SUM, LIS_MPI_COMM, ierr)
            call handle_mpi_error(ierr, &
                 'MPI_Allreduce call in USAF_superstatQC')
            superSigmaOSqr_pet(:) = 0

            superOErrScaleLength(:) = 0
            call MPI_Allreduce(superOErrScaleLength_pet, &
                 superOErrScaleLength, glbcr, &
                 MPI_REAL, MPI_SUM, LIS_MPI_COMM, ierr)
            call handle_mpi_error(ierr, &
                 'MPI_Allreduce call in USAF_superstatQC')
            superOErrScaleLength_pet(:) = 0
#endif

            ! Create the superobs
            do r = 1, LIS_rc%gnr(nest)
               do c = 1, LIS_rc%gnc(nest)
                  gindex = c + (r-1)*LIS_rc%gnc(nest)
                  if (superob_count(gindex) .gt. 2) then
                     superobs(gindex) = &
                          superobs(gindex) / superob_count(gindex)
                     superlat(gindex) = &
                          superlat(gindex) / superob_count(gindex)
                     superlon(gindex) = &
                          superlon(gindex) / superob_count(gindex)
                     superSigmaOSqr(gindex) = &
                          superSigmaOSqr(gindex) / superob_count(gindex)
                     superOErrScaleLength(gindex) = &
                          superOErrScaleLength(gindex) / superob_count(gindex)
                  else ! Superob either not possible or not needed
                     superobs(gindex) = MISSING
                  end if
               end do ! c
            end do ! r

         end if ! ipass .eq. 2

         ! Collect the final list of observations we will either reject or
         ! merge.
         if (ipass .eq. 3) then

#if (defined SPMD)
            call MPI_Barrier(LIS_MPI_COMM, ierr)
            call handle_mpi_error(ierr, &
                 'MPI_Barrier call in USAF_superstatQC')
            actions(:) = 0
            call MPI_Allreduce(actions_pet, actions, nobs,&
                 MPI_INTEGER, MPI_SUM, LIS_MPI_COMM, ierr)
            call handle_mpi_error(ierr, &
                 'MPI_Allreduce call in USAF_superstatQC')
            actions_pet(:) = 0
#endif
         end if ! ipass .eq. 3

      end do ! ipass

      ! Now update the observation structure with the rejected obs
      ! and the merged obs.  All processors do this to ensure
      ! identical modified structures.
      num_merged_obs = 0
      num_rejected_obs = 0
      do j = 1, nobs
         if (actions(j) .eq. -1) then
            this%qc(j) = QC_REJECT
            if (.not. silent_rejects_local) then
               write(LIS_logunit,*) &
                    '[INFO] superstatQC rejection j: ',j, &
                    ' net: ',trim(this%net(j)), &
                    ' platform: ',trim(this%platform(j)), &
                    ' lat: ',this%lat(j), &
                    ' lon: ',this%lon(j), &
                    ' obs: ',this%obs(j), &
                    ' back: ',this%back(j)
            endif
            num_rejected_obs = num_rejected_obs + 1
         else if (actions(j) .eq. 1) then
            this%qc(j) = QC_REJECT ! Was merged into superob
            num_merged_obs = num_merged_obs + 1
         end if
      end do ! j

      write(LIS_logunit,*) &
           '[INFO] superstatQC rejected ',num_rejected_obs, &
           ' obs and merged ',num_merged_obs,' obs'

      ! Finally, add the superobs to the data structure.
      num_superobs = 0
      do r = 1, LIS_rc%gnr(nest)
         do c = 1, LIS_rc%gnc(nest)
            gindex = c + (r-1)*LIS_rc%gnc(nest)
            if (superobs(gindex) .eq. MISSING) cycle

            call USAF_assignObsData(this,net_new,platform_new, &
                 superobs(gindex), &
                 superlat(gindex), &
                 superlon(gindex), &
                 superSigmaOSqr(gindex), &
                 superOErrScaleLength(gindex))

            num_superobs = num_superobs + 1
         end do ! c
      end do ! r

      write(LIS_logunit,*) &
           '[INFO] superstatQC created ',num_superobs,' super obs'

      ! Clean up
      deallocate(means)
      deallocate(actions)
      deallocate(superobs)
      deallocate(superlat)
      deallocate(superlon)
      deallocate(superSigmaOSqr)
      deallocate(superOErrScaleLength)
      deallocate(superob_count)

      deallocate(actions_pet)
      deallocate(superobs_pet)
      deallocate(superlat_pet)
      deallocate(superlon_pet)
      deallocate(superSigmaOSqr_pet)
      deallocate(superOErrScaleLength_pet)
      deallocate(superob_count_pet)

#if (defined SPMD)
      call MPI_Barrier(LIS_MPI_COMM, ierr)
      call handle_mpi_error(ierr, &
           'MPI_Barrier call in USAF_superstatQC')
      t2 = MPI_Wtime()
      write(LIS_logunit,*) &
           '[INFO] Elapsed time in superstatQC is ',t2 - t1,' seconds'
#endif

   end subroutine USAF_superstatQC

   !---------------------------------------------------------------------------
   ! QC checks for duplicate gage reports.  Based on Mahfouf et al (2007).
   !
   ! If duplicates are found for a particular station but all are identical,
   ! only one report is preserved and the rest are rejected.  Otherwise,
   ! if two different reports from the same station are found, a superob will
   ! be created if the spread is smaller than the observation error variance
   ! and both original reports will be rejected.  If more than two
   ! unique reports exist for the same station, all will be rejected.
   !
   ! EXCEPTION:  If the observation location is inconsistent, just save the
   ! first report.  This can happen if duplicate reports are received in
   ! different formats (e.g., BUFR and text) which have different precisions.
   ! Also, some observations are MOBL (the reporter literally is moving).
   subroutine USAF_dupQC(this)

      ! Imports
      use LIS_logMod, only: LIS_logunit
      use LIS_mpiMod

      ! Defaults
      implicit none

      ! Arguments
      type(USAF_ObsData), intent(inout) :: this

      ! Local variables
      type close_obs ! Linked list type
         integer :: ob_index
         type(close_obs), pointer :: next
      end type close_obs
      type(close_obs), pointer :: head, tail, new, ptr
      integer :: count_dups
      integer :: total_reject_count, total_merge_count, total_create_count
      real :: mean,back,newlat,newlon,sigmaOSqr,oErrScaleLength
      character(len=32) :: net
      character(len=32) :: platform
      real :: diff
      integer :: r,c,i
      integer :: nobs
      logical :: reject_all
      double precision :: t1, t2
      integer :: ierr
      logical :: location_issue

      nobs = this%nobs
      if (nobs .eq. 0) then
         write(LIS_logunit,*)&
              '[INFO] dupQC found no observations to test'
         return
      endif

#if (defined SPMD)
      call MPI_Barrier(LIS_MPI_COMM, ierr)
      call handle_mpi_error(ierr, &
           'MPI_Barrier call in USAF_dupQC')
      t1 = MPI_Wtime()
#endif

      total_reject_count = 0
      total_merge_count = 0
      total_create_count = 0

      do r = 1,nobs

         if (this%qc(r) .eq. QC_REJECT) cycle
         if (trim(this%net(r)) .eq. "SSMI") cycle
         if (trim(this%net(r)) .eq. "GEOPRECIP") cycle
         if (trim(this%net(r)) .eq. "CMORPH") cycle
         if (trim(this%net(r)) .eq. "IMERG") cycle

         ! Some CDMS obs are missing station IDs.  We will skip these
         if ( trim(this%net(r)) .eq. "CDMS" .and. &
              trim(this%platform(r)) .eq. "00000000") cycle

         ! Get count of duplicates of ob r
         count_dups = 0
         mean = 0
         nullify(head,tail) ! Initialize linked list

         do c = r+1,nobs
            if (this%qc(c) .eq. QC_REJECT) cycle
            if (trim(this%net(c)) .eq. "SSMI") cycle
            if (trim(this%net(c)) .eq. "GEOPRECIP") cycle
            if (trim(this%net(c)) .eq. "CMORPH") cycle
            if (trim(this%net(c)) .eq. "IMERG") cycle

            if ( trim(this%net(c)) .ne. trim(this%net(r))) cycle
            if ( trim(this%platform(c)) .ne. trim(this%platform(r))) cycle

            ! Some CDMS obs are missing station IDs.  We will skip these
            if ( trim(this%net(c)) .eq. "CDMS" .and. &
                 trim(this%platform(c)) .eq. "00000000") cycle

            ! Duplicate found.  Store in linked list.
            count_dups = count_dups + 1

            ! Add to linked list
            allocate(new)
            new%ob_index = c
            nullify(new%next)
            if (associated(head)) then
               tail%next => new
               tail => new
            else ! First entry
               head => new
               tail => new
            end if
         end do ! c

         ! If we have no duplicates, just move on.
         if (count_dups .eq. 0) cycle

         ! Handle the special case of inconsistent locations.  In this event,
         ! just save the first ob.
         location_issue = .false.
         if (count_dups .gt. 0) then

            ! First pass
            ptr => head
            do i = 1, count_dups
               diff = this%lat(ptr%ob_index) - this%lat(r)
               if (.not. diff > 0 .and. &
                    .not. diff < 0) location_issue = .true.
               diff = this%lon(ptr%ob_index) - this%lon(r)
               if (.not. diff > 0 .and. &
                    .not. diff < 0) location_issue = .true.
               if (location_issue) exit ! Get out of loop
               ptr => ptr%next
            end do ! i

            if (location_issue) then
               ! Toss all duplicate reports, just save the first one.
               ptr => head
               do i = 1, count_dups
                  this%qc(ptr%ob_index) = QC_REJECT
                  total_reject_count = total_reject_count + 1
                  write(LIS_logunit,*) &
                       '[INFO] dupQC rejecting ob w/ changing lat/lon  r: ', &
                       ptr%ob_index, &
                       ' net: ',trim(this%net(ptr%ob_index)), &
                       ' platform: ',trim(this%platform(ptr%ob_index)), &
                       ' lat: ',this%lat(ptr%ob_index), &
                       ' lon: ',this%lon(ptr%ob_index), &
                       ' obs: ',this%obs(ptr%ob_index), &
                       ' back: ',this%back(ptr%ob_index)
                  write(LIS_logunit,*) &
                       '------------------------------------------------------'
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
               diff = this%obs(ptr%ob_index) - this%obs(r)
               if (.not. diff > 0 .and. &
                    .not. diff < 0) then
                  this%qc(ptr%ob_index) = QC_REJECT
                  total_reject_count = total_reject_count + 1
                  write(LIS_logunit,*) &
                       '[INFO] dupQC rejecting exact duplicate1 ob r: ', &
                       ptr%ob_index, &
                       ' net: ',trim(this%net(ptr%ob_index)), &
                       ' platform: ',trim(this%platform(ptr%ob_index)), &
                       ' lat: ',this%lat(ptr%ob_index), &
                       ' lon: ',this%lon(ptr%ob_index), &
                       ' obs: ',this%obs(ptr%ob_index), &
                       ' back: ',this%back(ptr%ob_index)
                  write(LIS_logunit,*) &
                       '------------------------------------------------------'

               else
                  reject_all = .true.
                  exit ! out of i loop
               end if
               ptr => ptr%next
            end do ! i

            if (reject_all) then
               this%qc(r) = QC_REJECT
               total_reject_count = total_reject_count + 1
               write(LIS_logunit,*) &
                    '[INFO] dupQC rejecting ob1 r: ', &
                    r, &
                    ' net: ',trim(this%net(r)), &
                    ' platform: ',trim(this%platform(r)), &
                    ' lat: ',this%lat(r), &
                    ' lon: ',this%lon(r), &
                    ' obs: ',this%obs(r), &
                    ' back: ',this%back(r)

               ptr => head
               do i = 1, count_dups
                  this%qc(ptr%ob_index) = QC_REJECT
                  total_reject_count = total_reject_count + 1
                  write(LIS_logunit,*) &
                       '[INFO] dupQC rejecting ob1 r: ', &
                       ptr%ob_index, &
                       ' net: ',trim(this%net(ptr%ob_index)), &
                       ' platform: ',trim(this%platform(ptr%ob_index)), &
                       ' lat: ',this%lat(ptr%ob_index), &
                       ' lon: ',this%lon(ptr%ob_index), &
                       ' obs: ',this%obs(ptr%ob_index), &
                       ' back: ',this%back(ptr%ob_index)

                  ptr => ptr%next
               end do ! i
               write(LIS_logunit,*) &
                    '------------------------------------------------------'
            end if ! reject_all
         end if ! count_dups .gt. 1 .and. .not. location_issue

         ! If we have exactly one duplicate: reject duplicate
         ! if it is an exact copy; otherwise, attempt superob.
         if (count_dups .eq. 1 .and. .not. location_issue) then
            ptr => head
            diff = this%obs(ptr%ob_index) - this%obs(r)
            if (.not. diff < 0 .and. &
                 .not. diff > 0) then
               this%qc(ptr%ob_index) = QC_REJECT
               total_reject_count = total_reject_count + 1
               write(LIS_logunit,*) &
                    '[INFO] dupQC rejecting exact duplicate2 ob r: ', &
                    ptr%ob_index, &
                    ' net: ',trim(this%net(ptr%ob_index)), &
                    ' platform: ',trim(this%platform(ptr%ob_index)), &
                    ' lat: ',this%lat(ptr%ob_index), &
                    ' lon: ',this%lon(ptr%ob_index), &
                    ' obs: ',this%obs(ptr%ob_index), &
                    ' back: ',this%back(ptr%ob_index)

               write(LIS_logunit,*) &
                    '------------------------------------------------------'

            else if (diff*diff .gt. this%sigmaOSqr(r)) then
               this%qc(r) = QC_REJECT
               total_reject_count = total_reject_count + 1
               write(LIS_logunit,*) &
                    '[INFO] dupQC rejecting2 ob r: ', &
                    r, &
                    ' net: ',trim(this%net(r)), &
                    ' platform: ',trim(this%platform(r)), &
                    ' lat: ',this%lat(r), &
                    ' lon: ',this%lon(r), &
                    ' obs: ',this%obs(r), &
                    ' back: ',this%back(r)

               this%qc(ptr%ob_index) = QC_REJECT
               total_reject_count = total_reject_count + 1

               write(LIS_logunit,*) &
                    '[INFO] dupQC rejecting2 ob r: ', &
                    ptr%ob_index, &
                    ' net: ',trim(this%net(ptr%ob_index)), &
                    ' platform: ',trim(this%platform(ptr%ob_index)), &
                    ' lat: ',this%lat(ptr%ob_index), &
                    ' lon: ',this%lon(ptr%ob_index), &
                    ' obs: ',this%obs(ptr%ob_index), &
                    ' back: ',this%back(ptr%ob_index)

               write(LIS_logunit,*) &
                    '------------------------------------------------------'

            else
               mean = 0.5 * (this%obs(ptr%ob_index) + this%obs(r))

               write(LIS_logunit,*) &
                    '[INFO] dupQC will create superob from r: ', &
                    r, &
                    ' net: ',trim(this%net(r)), &
                    ' platform: ',trim(this%platform(r)), &
                    ' lat: ',this%lat(r), &
                    ' lon: ',this%lon(r), &
                    ' obs: ',this%obs(r), &
                    ' back: ',this%back(r)

               write(LIS_logunit,*) &
                    '[INFO] dupQC will create superob from r: ', &
                    ptr%ob_index, &
                    ' net: ',trim(this%net(ptr%ob_index)), &
                    ' platform: ',trim(this%platform(ptr%ob_index)), &
                    ' lat: ',this%lat(ptr%ob_index), &
                    ' lon: ',this%lon(ptr%ob_index), &
                    ' obs: ',this%obs(ptr%ob_index), &
                    ' back: ',this%back(ptr%ob_index)

               mean = mean
               back = this%back(r)
               newlat = this%lat(r)
               newlon = this%lon(r)
               sigmaOSqr = this%sigmaOSqr(r)
               oErrScaleLength = this%oErrScaleLength(r)

               ! EMK Bug fix:  Make sure net and platform are set for the
               ! superob.  Print these values.
               net = trim(this%net(ptr%ob_index))
               platform = trim(this%platform(ptr%ob_index))
               write(LIS_logunit,*) &
                    '[INFO] dupQC new superob is: ', &
                    ' net: ',trim(net), &
                    ' platform: ',trim(platform), &
                    ' obs: ',mean
               write(LIS_logunit,*) &
                    '------------------------------------------------------'

               call USAF_assignObsData(this,net,platform,mean,newlat,newlon,&
                    sigmaOSqr,oErrScaleLength,back=back)

               total_create_count = total_create_count + 1

               ! Reject originals
               this%qc(r) = QC_REJECT
               this%qc(ptr%ob_index) = QC_REJECT
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

      end do ! r

      write(LIS_logunit,*) &
           '[INFO] dupQC rejected ',total_reject_count,' obs and merged ', &
           total_merge_count
      write(LIS_logunit,*) &
           '[INFO] dupQC created ',total_create_count,' super obs'

#if (defined SPMD)
      call MPI_Barrier(LIS_MPI_COMM, ierr)
      call handle_mpi_error(ierr, &
           'MPI_Barrier call in USAF_dupQC')
      t2 = MPI_Wtime()
      write(LIS_logunit,*) &
           '[INFO] Elapsed time in dupQC is ',t2 - t1,' seconds'
#endif

   end subroutine USAF_dupQC

   !---------------------------------------------------------------------------
   ! Selects new ESMF PET ID based on current value and increment.
   ! If resulting ID is out of bounds, it is reset to the nearest bound and
   ! the increment is reversed.  Resulting behavior is a 1-D "ping-pong"
   ! of values between 0 and (LIS_npes-1).
   subroutine update_pet(pet,pet_incr)
      use LIS_coreMod, only: LIS_npes
      implicit none
      integer,intent(inout) :: pet
      integer,intent(inout) :: pet_incr
      pet = pet + pet_incr
      if (pet .ge. LIS_npes) then
         pet = LIS_npes - 1
         pet_incr = -1
      else if (pet .lt. 0) then
         pet = 0
         pet_incr = 1
      end if
   end subroutine update_pet

   !---------------------------------------------------------------------------
   ! Changes negative values in array to zero, and summarizes number of
   ! values so changed.
   subroutine reset_negative_values(nest,mrgp)

      ! Imports
      use LIS_coreMod, only: LIS_rc
      use LIS_logMod, only: LIS_logunit

      ! Arguments
      integer, intent(in) :: nest
      real, intent(inout) :: mrgp(LIS_rc%lnc(nest),LIS_rc%lnr(nest))

      ! Local variables
      integer :: ifix
      integer :: c,r

      ifix = 0
      do r = 1, LIS_rc%lnr(nest)
         do c = 1, LIS_rc%lnc(nest)
            if (mrgp(c,r) .lt. 0) then
               mrgp(c,r) = 0
               ifix = ifix + 1
            end if
         end do ! c
      end do ! r
      if (ifix > 0) then
         write(LIS_logunit,6000) ifix
6000     format (/, 1x, 55('-'), &
              /, 3x, '[INFO] routine reset_negative_values:',&
              /, 5x, '# of pts to which negative values were set to zero = ', &
              i6, /, 1x, 55('-'))
      end if

   end subroutine reset_negative_values

   !---------------------------------------------------------------------------
   ! Reject obs that differ "too much" from background field.  Threshold based
   ! on sum of observation and background error variances.  From Lopez (2013).
   ! This assumes both background and observations are unbiased, and large
   ! difference implies gross error in observation.
   subroutine USAF_backQC(this,sigmaBSqr,silent_rejects)

      ! Imports
      use LIS_logMod, only: LIS_logunit, LIS_endrun
      use LIS_mpiMod

      ! Defaults
      implicit none

      ! Arguments
      type(USAF_ObsData), intent(inout) :: this
      real, intent(in) :: sigmaBSqr
      logical,optional,intent(in) :: silent_rejects

      ! Local variables
      integer :: nobs
      real :: errorThresh
      real :: absDiff
      integer :: r
      integer :: reject_count
      integer :: ierr
      double precision :: t1, t2
      logical :: silent_rejects_local

      nobs = this%nobs
      if (nobs .eq. 0) then
         write(LIS_logunit,*)&
              '[INFO] backQC found no observations to test'
         return
      endif

      silent_rejects_local = .false.
      if (present(silent_rejects)) then
         silent_rejects_local = silent_rejects
      end if

#if (defined SPMD)
      call MPI_Barrier(LIS_MPI_COMM, ierr)
      call handle_mpi_error(ierr, &
           'MPI_Barrier call in USAF_backQC')
      t1 = MPI_Wtime()
#endif

      !write(LIS_logunit,*)'EMK: nobs, sigmaBSqr = ', nobs, sigmaBSqr
      !flush(LIS_Logunit)

      reject_count = 0
      do r = 1,nobs

         ! Skip bad data
         if ( this%qc(r) .eq. QC_REJECT) cycle

         errorThresh = 4*sqrt(sigmaBSqr + this%sigmaOSqr(r))
         absDiff = abs(this%obs(r) - this%back(r))

         if (absDiff .gt. errorThresh) then
            this%qc(r) = QC_REJECT

            reject_count = reject_count + 1

            if (.not. silent_rejects_local) then
               write(LIS_logunit,*) &
                    '[INFO] backQC rejecting observation i: ',r, &
                    ' net: ',trim(this%net(r)), &
                    ' platform: ',trim(this%platform(r)), &
                    ' lat: ',this%lat(r), &
                    ' lon: ',this%lon(r), &
                    ' obs: ',this%obs(r), &
                    ' back: ',this%back(r), &
                    ' abs diff: ', abs(this%obs(r) - this%back(r)), &
                    ' errorThresh ', errorThresh
            end if
         end if

      end do ! r

      write(LIS_logunit,*)&
           '[INFO] backQC rejected ',reject_count,' observations'

#if (defined SPMD)
      call MPI_Barrier(LIS_MPI_COMM, ierr)
      call handle_mpi_error(ierr, &
           'MPI_Barrier call in USAF_backQC')
      t2 = MPI_Wtime()
      write(LIS_logunit,*) &
           '[INFO] Elapsed time in backQC is ',t2 - t1,' seconds'
#endif

   end subroutine USAF_backQC

   !---------------------------------------------------------------------------
   ! Checks if observation network is recognized as a gauge.
   logical function is_gauge(net)
      implicit none
      character(len=32), intent(in) :: net
      logical :: answer
      answer = .false.
      if (trim(net) .eq. "AMIL") answer = .true.
      if (trim(net) .eq. "CANA") answer = .true.
      if (trim(net) .eq. "FAA") answer = .true.
      if (trim(net) .eq. "ICAO") answer = .true.
      if (trim(net) .eq. "WMO") answer = .true.
      if (trim(net) .eq. "MOBL") answer = .true.
      if (trim(net) .eq. "SUPERGAGE") answer = .true.
      ! Handle reformatted CDMS data that are missing the network type.
      if (trim(net) .eq. "CDMS") answer = .true.
      is_gauge = answer
   end function is_gauge

   !---------------------------------------------------------------------------
   ! Dummy function for establishing a surface station is uncorrelated.
   ! Simple returns "true".  This is passed to some of the generic
   ! Bratseth routines that need to know which reports in a collection
   ! have correlated errors.  When analyzing screen-level variables with
   ! surface stations, all observations should have uncorrelated errors.
   logical function is_stn(net)
      implicit none
      character(len=32), intent(in) :: net
      logical :: answer
      answer = .true.
      is_stn = answer
   end function is_stn

   !---------------------------------------------------------------------------
   ! Checks if observation "network" is recognized as SSMI retrievals.
   logical function is_ssmi(net)
      implicit none
      character(len=32), intent(in) :: net
      logical :: answer
      answer = .false.
      if (trim(net) .eq. "SSMI") answer = .true.
      is_ssmi = answer
   end function is_ssmi

   !---------------------------------------------------------------------------
   ! Checks if observation "network" is recognized as GEOPRECIP retrievals.
   logical function is_geoprecip(net)
      implicit none
      character(len=32), intent(in) :: net
      logical :: answer
      answer = .false.
      if (trim(net) .eq. "GEOPRECIP") answer = .true.
      is_geoprecip = answer
   end function is_geoprecip

   !---------------------------------------------------------------------------
   ! Checks if observation "network" is recognized as CMORPH estimate.
   logical function is_cmorph(net)
      implicit none
      character(len=32), intent(in) :: net
      logical :: answer
      answer = .false.
      if (trim(net) .eq. "CMORPH") answer = .true.
      is_cmorph = answer
   end function is_cmorph

   !---------------------------------------------------------------------------
   ! Checks if observation "network" is recognized as IMERG retrievals.
   logical function is_imerg(net)
      implicit none
      character(len=32), intent(in) :: net
      logical :: answer
      answer = .false.
      if (trim(net) .eq. "IMERG") answer = .true.
      is_imerg = answer
   end function is_imerg

   !---------------------------------------------------------------------------
   ! Initializes 2d hash table for storing observations (actually array index
   ! values) based on LIS grid box row and column.
   subroutine init_hash2d(imax,jmax,hash2d)

      ! Defaults
      implicit none

      ! Arguments
      integer, intent(in) :: imax
      integer, intent(in) :: jmax
      type(hash_list), allocatable, intent(out) :: hash2d(:,:)

      ! Local variables
      integer :: i,j

      allocate(hash2d(imax,jmax))

      do j = 1, jmax
         do i = 1, imax
            hash2d(i,j)%obindex = MISSING
            nullify(hash2d(i,j)%next)
         end do ! i
      end do ! j

   end subroutine init_hash2d

   !---------------------------------------------------------------------------
   ! Inserts array index for an observation into the 2d hash table organizing
   ! these indices by LIS grid row and column.
   subroutine insert_hash2d(imax,jmax,hash2d,i,j,obindex)

      ! Imports

      ! Defaults
      implicit none

      ! Arguments
      integer, intent(in) :: imax
      integer, intent(in) :: jmax
      type(hash_list), target, intent(inout) :: hash2d(imax,jmax)
      integer, intent(in) :: i
      integer, intent(in) :: j
      integer, intent(in) :: obindex

      ! Local variables
      type(hash_list), pointer :: new_node
      type(hash_list), pointer :: node

      nullify(node, new_node)

      node => hash2d(i,j)
      do
         if (node%obindex .eq. MISSING) then
            node%obindex = obindex
            exit
         else
            if (.not. associated(node%next)) then
               allocate(new_node)
               new_node%obindex = obindex
               nullify(new_node%next)
               node%next => new_node
               exit
            else
               node => node%next
            end if
         end if
      end do

   end subroutine insert_hash2d

   !---------------------------------------------------------------------------
   subroutine insert_hash2d_array(imax,jmax,hash2d,nobs,cols,rows)

      ! Imports

      ! Defaults
      implicit none

      ! Arguments
      integer, intent(in) :: imax
      integer, intent(in) :: jmax
      type(hash_list), target, intent(inout) :: hash2d(imax,jmax)
      integer, intent(in) :: nobs
      integer, intent(in) :: cols(nobs)
      integer, intent(in) :: rows(nobs)

      ! Local variables
      integer :: j
      type(hash_list), pointer :: new_node
      type(hash_list), pointer :: node

      do j = 1, nobs
         if (cols(j) .eq. 0) cycle
         if (rows(j) .eq. 0) cycle
         nullify(node, new_node)
         node => hash2d(cols(j),rows(j))
         do
            if (node%obindex .eq. MISSING) then
               node%obindex = j
               exit
            else
               if (.not. associated(node%next)) then
                  allocate(new_node)
                  new_node%obindex = j
                  nullify(new_node%next)
                  node%next => new_node
                  exit
               else
                  node => node%next
               end if
            end if
         end do
      end do
   end subroutine insert_hash2d_array

   !---------------------------------------------------------------------------
   ! Return list of array indices for all observations within requested
   ! LIS grid row and column.  Leverages 2d hash table.
   subroutine get_hash2d(imax,jmax,hash2d,i,j,nobs_local,obindex_vector)

      ! Defaults
      implicit none

      ! Arguments
      integer, intent(in) :: imax
      integer, intent(in) :: jmax
      type(hash_list), target, intent(in) :: hash2d(imax,jmax)
      integer, intent(in) :: i
      integer, intent(in) :: j
      integer, intent(out) :: nobs_local
      integer, allocatable, intent(out) :: obindex_vector(:)

      ! Local variables
      type(hash_list), pointer :: node
      integer :: k

      nullify(node)

      ! First, find out how many obs are in this grid box
      nobs_local = 0
      if (hash2d(i,j)%obindex .eq. MISSING) return

      node => hash2d(i,j)
      nobs_local = nobs_local + 1
      do
         if (.not. associated(node%next)) then
            exit
         else
            node => node%next
            nobs_local = nobs_local + 1
         end if
      end do

      allocate(obindex_vector(nobs_local))

      ! Second pass:  Get all the ob gindex values
      node => hash2d(i,j)
      k = 1
      obindex_vector(k) = node%obindex
      do
         if (.not. associated(node%next)) then
            exit
         else
            node => node%next
            k = k + 1
            obindex_vector(k) = node%obindex
         end if
      end do

   end subroutine get_hash2d

   !---------------------------------------------------------------------------
   ! Dismantles 2d hash table of array indices of observations in each
   ! LIS grid box row and column.  The linked lists in the hash table are
   ! deallocated, but hash2d must be deallocated separately after calling
   ! this subroutine.
   subroutine destroy_hash2d(imax,jmax,hash2d)

      ! Defaults
      implicit none

      ! Arguments
      integer, intent(in) :: imax
      integer, intent(in) :: jmax
      type(hash_list), target, intent(inout) :: hash2d(imax,jmax)

      ! Local variables
      type(hash_list), pointer :: node, next, first
      integer :: i,j

      nullify(node, next, first)

      do j = 1, jmax
         do i = 1, imax
            ! Very first node in list must be preserved until we deallocate
            ! hash2d. But all nodes beyond that can be deallocated one by one.
            first => hash2d(i,j)
            first%obindex = MISSING
            if (.not. associated(first%next)) cycle
            node => first%next
            nullify(first%next)
            do
               if (associated(node%next)) then
                  next => node%next
                  deallocate(node)
                  node => next
               else ! Last node in list
                  deallocate(node)
                  nullify(node)
                  nullify(next)
                  exit
               end if
            end do
         end do ! i
      end do ! j

   end subroutine destroy_hash2d

   !---------------------------------------------------------------------------
   ! Locates the "good" observations in the ObsData structure, and populates
   ! a 2d hash table with the array indices based on LIS grid row and column.
   subroutine build_hash2d(this,nest,imax,jmax,hash2d)

      ! Imports
      use LIS_coreMod, only: LIS_rc
      use LIS_mpiMod

      ! Defaults
      implicit none

      ! Arguments
      type(USAF_ObsData), intent(inout) :: this
      integer,intent(in) :: nest
      integer,intent(out) :: imax
      integer,intent(out) :: jmax
      type(hash_list), allocatable, intent(out) :: hash2d(:,:)

      ! Local variables
      integer :: nobs
      integer, allocatable :: cols(:), rows(:)

      nobs = this%nobs
      if (nobs .eq. 0) return

      ! Get the columns and rows of the observations in the LIS grid
      call find_LIS_cols_rows(this, nest, cols, rows)

      ! Loop through the collected col, row information and add to the
      ! hash2d table.
      call init_hash2d(LIS_rc%gnc(nest), LIS_rc%gnr(nest), hash2d)
      imax = LIS_rc%gnc(nest) ! Output
      jmax = LIS_rc%gnr(nest) ! Output
      call insert_hash2d_array(imax, jmax, hash2d, nobs, cols, rows)
!       do j = 1,nobs
!          c = cols(j)
!          if (c .eq. 0) cycle
!          r = rows(j)
!          if (r .eq. 0) cycle
!          call insert_hash2d(LIS_rc%gnc(nest),LIS_rc%gnr(nest),hash2d,c,r,j)
!       end do

      ! Clean up
      if (allocated(rows)) deallocate(rows)
      if (allocated(cols)) deallocate(cols)

   end subroutine build_hash2d

   !---------------------------------------------------------------------------
   ! Find neighbors around a designated grid point close enough to have
   ! partial correlation in the background field.  Assumes lat/lon grid,
   ! and leverages a 2d hash table filled with observation array indices.
   !
   ! The neighbors are indicated by one set of lower and upper rows, and
   ! possibly two sets of left and right columns.  The optional second set
   ! is necessary if a global lat/lon LIS domain is used, and the designated
   ! grid point is near the eastern or western lateral boundary; in this case,
   ! neighbors may "wrap around" the LIS domain and be located near the
   ! opposite lateral boundary.
   !
   ! "Neighbors" are determined by comparing great circle distance with a
   ! radius of influence.
   subroutine find_gridpt_neighbors(c,r,nest,max_dist, &
        lr,ur,lc1,rc1,lc2,rc2)

      ! Imports
      use LIS_coreMod, only:  LIS_domain, LIS_rc

      ! Defaults
      implicit none

      ! Arguments
      integer, intent(in) :: c
      integer, intent(in) :: r
      integer, intent(in) :: nest
      real, intent(in) :: max_dist
      integer, intent(out) :: lr
      integer, intent(out) :: ur
      integer, intent(out) :: lc1
      integer, intent(out) :: rc1
      integer, intent(out) :: lc2
      integer, intent(out) :: rc2

      ! Local variables
      integer :: gindex
      real :: ctrlat, tmplat
      real :: ctrlat_lower,ctrlon_lower
      real :: ctrlat_upper,ctrlon_upper
      integer :: lc1_upper,rc1_upper,lc2_upper,rc2_upper
      integer :: lc1_lower,rc1_lower,lc2_lower,rc2_lower
      real :: dist

      ! Find latitude bounds.  Assumes lat/lon domain
      gindex = 1 + (r-1)*LIS_rc%gnc(nest)
      ctrlat = LIS_domain(nest)%glat(gindex)

      ! First, lower bound
      do lr = r, 1, -1
         gindex = 1 + (lr-1)*LIS_rc%gnc(nest)
         tmplat = LIS_domain(nest)%glat(gindex)
         dist = great_circle_distance(ctrlat, 0., tmplat, 0.)
         if (dist .gt. max_dist) exit
      end do ! lr
      lr = max(lr-1,1)

      ! Next, upper bound
      do ur = r, LIS_rc%gnr(nest), 1
         gindex = 1 + (ur-1)*LIS_rc%gnc(nest)
         tmplat = LIS_domain(nest)%glat(gindex)
         dist = great_circle_distance(ctrlat, 0., tmplat, 0.)
         if (dist .gt. max_dist) exit
      end do ! ur
      ur = min(ur+1,LIS_rc%gnr(nest))

      ! Find longitude bounds.  Assumes lat/lon domain
      ! First, use lower latitude.
      gindex = c + (lr-1)*LIS_rc%gnc(nest)
      ctrlat_lower = LIS_domain(nest)%glat(gindex)
      ctrlon_lower = LIS_domain(nest)%glon(gindex)
      call find_gridpt_neighbors_leftlon(c,lr,nest,ctrlat_lower,ctrlon_lower, &
           max_dist,lc1_lower,rc2_lower)
      call find_gridpt_neighbors_rightlon(c,lr,nest,ctrlat_lower,ctrlon_lower,&
           max_dist,rc1_lower,lc2_lower)

      ! Next, use upper latitude.
      gindex = c + (ur-1)*LIS_rc%gnc(nest)
      ctrlat_upper = LIS_domain(nest)%glat(gindex)
      ctrlon_upper = LIS_domain(nest)%glon(gindex)
      call find_gridpt_neighbors_leftlon(c,ur,nest,ctrlat_upper,ctrlon_upper,&
           max_dist,lc1_upper,rc2_upper)
      call find_gridpt_neighbors_rightlon(c,ur,nest,ctrlat_lower,ctrlon_lower,&
           max_dist,rc1_upper,lc2_upper)

      ! Use the larger of the two bounds (associated with poleward latitude)
      if (abs(ctrlat_upper) .gt. abs(ctrlat_lower)) then
         lc1 = lc1_upper
         rc2 = rc2_upper
         rc1 = rc1_upper
         lc2 = lc2_upper
      else
         lc1 = lc1_lower
         rc2 = rc2_lower
         rc1 = rc1_lower
         lc2 = lc2_lower
      end if

      ! EMK TEST
      if (rc2 .eq. -1) then
         print*,'EMK ERROR: c,r,lr,ur,lc1,rc2,rc1,lc2 = ', &
              c,r,lr,ur,lc1,rc2,rc1,lc2
         print*,'EMK ERROR: lc1_upper,rc2_upper,rc1_upper,lc2_upper = ', &
              lc1_upper,rc2_upper,rc1_upper,lc2_upper
         print*,'EMK ERROR: lc1_lower,rc2_lower,rc1_lower,lc2_lower = ', &
              lc1_lower,rc2_lower,rc1_lower,lc2_lower
         stop 1
      end if

   end subroutine find_gridpt_neighbors

   !---------------------------------------------------------------------------
   ! This finds neighbors to the left (west) of the designated grid box.
   ! Accounts for "wrap-around" effect of global lat/lon grid by returning
   ! two index bounds, one truly left (west), and one possibly right (east).
   !
   ! "Neighbors" are determined by comparing great circle distance with a
   ! radius of influence.
   subroutine find_gridpt_neighbors_leftlon(c,r,nest,ctrlat,ctrlon,max_dist, &
        lc1,rc2)

      ! Imports
      use LIS_coreMod, only: LIS_rc, LIS_domain

      ! Defaults
      implicit none

      ! Arguments
      integer, intent(in) :: c
      integer, intent(in) :: r
      integer, intent(in) :: nest
      real, intent(in) :: ctrlat
      real, intent(in) :: ctrlon
      real, intent(in) :: max_dist
      integer, intent(out) :: lc1
      integer, intent(out) :: rc2

      ! Local variables
      logical :: found
      integer :: gindex
      real :: tmplon
      real :: dist

      found = .false.
      do lc1 = c, 1, -1
         gindex = lc1 + (r-1)*LIS_rc%gnc(nest)
         tmplon = LIS_domain(nest)%glon(gindex)
         dist = great_circle_distance(ctrlat, ctrlon, ctrlat, tmplon)
         if (dist .gt. max_dist) then
            found = .true.
            exit
         end if
      end do ! lc1
      if (.not. found) then
         do rc2 = LIS_rc%gnc(nest), c, -1
            gindex = rc2 + (r-1)*LIS_rc%gnc(nest)
            tmplon = LIS_domain(nest)%glon(gindex)
            dist = great_circle_distance(ctrlat, ctrlon, ctrlat, tmplon)
            if (dist .gt. max_dist) then
               found = .true.
               exit
            end if
         end do ! rc2
         if (found) then
            rc2 = rc2-1
            lc1 = 1
         else ! Entire latitudinal band is in range
            lc1 = 1
            rc2 = 0
         end if
      else
         lc1 = max(lc1-1,1)
         rc2 = 0
      end if

   end subroutine find_gridpt_neighbors_leftlon

   !---------------------------------------------------------------------------
   ! This finds neighbors to the right (east) of the designated grid box.
   ! Accounts for "wrap-around" effect of global lat/lon grid by returning
   ! two index bounds, one truly right (east), and one possibly left (west).
   !
   ! "Neighbors" are determined by comparing great circle distance with a
   ! radius of influence.
   subroutine find_gridpt_neighbors_rightlon(c,r,nest,ctrlat,ctrlon,max_dist, &
        rc1,lc2)

      ! Imports
      use LIS_coreMod, only: LIS_rc, LIS_domain

      ! Defaults
      implicit none

      ! Arguments
      integer, intent(in) :: c
      integer, intent(in) :: r
      integer, intent(in) :: nest
      real, intent(in) :: ctrlat
      real, intent(in) :: ctrlon
      real, intent(in) :: max_dist
      integer, intent(out) :: rc1
      integer, intent(out) :: lc2

      ! Local variables
      logical :: found
      integer :: gindex
      real :: tmplon
      real :: dist

      ! Find right longitude bound.  It may wrap around the grid, so
      ! two indices may be found.
      found = .false.
      do rc1 = c, LIS_rc%gnc(nest)
         gindex = rc1 + (r-1)*LIS_rc%gnc(nest)
         tmplon = LIS_domain(nest)%glon(gindex)
         dist = great_circle_distance(ctrlat, ctrlon, ctrlat, tmplon)
         if (dist .gt. max_dist) then
            found = .true.
            exit
         end if
      end do ! rc1
      if (.not. found) then
         do lc2 = 1, c
            gindex = lc2 + (r-1)*LIS_rc%gnc(nest)
            tmplon = LIS_domain(nest)%glon(gindex)
            dist = great_circle_distance(ctrlat, ctrlon, ctrlat, tmplon)
            if (dist .gt. max_dist) then
               found = .true.
               exit
            end if
         end do ! rc2
         if (found) then
            lc2 = lc2+1
            rc1 = LIS_rc%gnc(nest)
         else ! Entire latitudinal band is in range
            rc1 = LIS_rc%gnc(nest)
            lc2 = 0
         end if
      else
         rc1 = min(rc1+1,LIS_rc%gnc(nest))
         lc2 = 0
      end if

   end subroutine find_gridpt_neighbors_rightlon

   !---------------------------------------------------------------------------
   ! Collects array indices of all observations surrounding a given LIS
   ! grid row and column.  Uses 2d hash table storing these indices according
   ! to global grid position.  The top (north) and bottom (south) bounds of
   ! the search grid must be provided, plus two sets of left (west) and right
   ! (east) bounds to account for possible "wrap-around" the global lat/lon
   ! domain.
   subroutine get_neighbor_obs(cmax,rmax,hash2d,lr,ur,lc1,rc1,lc2,rc2,nest, &
        nobs_neighbors, iobs_neighbors_vector)

      ! Imports
      use LIS_coreMod, only: LIS_rc

      ! Defaults
      implicit none

      ! Arguments
      integer, intent(in) :: cmax
      integer, intent(in) :: rmax
      type(hash_list), intent(in) :: hash2d(cmax,rmax)
      integer, intent(in) :: lr
      integer, intent(in) :: ur
      integer, intent(in) :: lc1
      integer, intent(in) :: rc1
      integer, intent(in) :: lc2
      integer, intent(in) :: rc2
      integer, intent(in) :: nest
      integer, intent(out) :: nobs_neighbors
      integer, allocatable, intent(out) :: iobs_neighbors_vector(:)

      ! Local variables
      integer :: c,r,jpass,ipass
      integer :: start_c, end_c
      integer :: obcount_cr
      integer, allocatable :: iobs_cr_vector(:)
      integer :: icount
      integer :: i

      ! jpass=1 is for counting obs, so we can allocate the return vector.
      ! jpass=2 is for copying data to the return vector.
      do jpass = 1,2
         icount = 0
         do r = lr,ur
            ! ipass=1 is the normal search for neighbors, from left
            ! and right longitude bounds.
            ! ipass=2 and ipass=3 are special cases where the central grid
            ! point is near the meridian lateral boundary, and neighbors
            ! may exist on the other side of that boundary.
            do ipass = 1,3
               if (ipass .eq. 1) then
                  start_c = lc1
                  end_c = rc1
               else if (ipass .eq. 2) then
                  if (rc2 .eq. 0) cycle
                  start_c = rc2
                  end_c = LIS_rc%gnc(nest)
               else if (ipass .eq. 3) then
                  if (lc2 .eq. 0) cycle
                  start_c = 1
                  end_c = lc2
               end if

!                !EMK DEBUG
!                if (start_c < 1) then
!                   print*,'EMK ERROR: jpass,ipass,lc1,rc2,start_c = ', &
!                        jpass,ipass,lc1,rc2,start_c
!                   stop 1
!                end if

               do c = start_c, end_c
                  call get_hash2d(LIS_rc%gnc(nest),LIS_rc%gnr(nest),&
                       hash2d,c,r,&
                       obcount_cr,iobs_cr_vector)

                  if (jpass .eq. 2 .and. obcount_cr .gt. 0) then
                     do i = 1, obcount_cr
                        iobs_neighbors_vector(icount+i) = iobs_cr_vector(i)
                     end do ! i
                  end if
                  icount = icount + obcount_cr

                  if (allocated(iobs_cr_vector)) deallocate(iobs_cr_vector)

               end do ! c
            end do ! ipass
         end do ! r

         if (jpass .eq. 1) then
            nobs_neighbors = icount
            if (nobs_neighbors .eq. 0) exit
            allocate(iobs_neighbors_vector(nobs_neighbors))
            icount = 0
         end if

      end do ! jpass

   end subroutine get_neighbor_obs


   !---------------------------------------------------------------------------
   ! Find the column and rows of "good" observations in LIS domain.
   subroutine find_LIS_cols_rows(this,nest,cols,rows)

      ! Imports
      use LIS_coreMod, only: LIS_rc, LIS_domain, LIS_localPet
      use LIS_logMod, only: LIS_logunit
      use LIS_mpiMod

      ! Defaults
      implicit none

      ! Arguments
      type(USAF_ObsData), intent(inout) :: this
      integer,intent(in) :: nest
      integer, allocatable,intent(out) :: cols(:), rows(:)

      ! Local variables
      real :: dlat, dlon, ctrlat, ctrlon
      integer :: nobs
      integer, allocatable :: cols_pet(:), rows_pet(:)
      integer :: pet, pet_incr
      logical :: found
      integer :: j,r,c,gindex
      integer :: ierr

      nobs = this%nobs
      if (nobs .eq. 0) return

      dlat = LIS_domain(nest)%lisproj%dlat
      dlon = LIS_domain(nest)%lisproj%dlon

      ! Prepare to look through all observations and find the LIS grid
      ! box they belong in based on lat/lon.  This work is parallelized
      ! due to the shear size of some satellite data (e.g., IMERG)
      allocate(cols_pet(nobs))
      cols_pet = 0
      allocate(rows_pet(nobs))
      rows_pet = 0
      pet = -1
      pet_incr = 1

#if (defined SPMD)
      call MPI_Barrier(LIS_MPI_COMM, ierr)
      call handle_mpi_error(ierr, &
           'MPI_Barrier call in find_LIS_cols_rows')
#endif
      do j = 1, nobs

         if (this%qc(j) .eq. QC_REJECT) cycle

         ! See which PET is responsible for this ob
         call update_pet(pet,pet_incr)
         if (pet .ne. LIS_localPet) cycle

         ! First, latitude
         found = .false.
         do r = 1, LIS_rc%gnr(nest)
            gindex = 1 + (r-1)*LIS_rc%gnc(nest)
            ctrlat = LIS_domain(nest)%glat(gindex)

            ! EMK 8 June 2017...Imprecise floating point arithmetic
            ! can cause an ob with a latitude exactly lying on a grid box
            ! edge to "fall between the cracks".  To avoid this, we generally
            ! only check if the ob is south of the northern grid box edge.
            ! An exception is along the first row, where we need to see
            ! if the ob is south of the LIS domain.
            if (r .eq. 1) then
               if (this%lat(j) .lt. (ctrlat - (0.5*dlat))) cycle
            end if
            if (this%lat(j) .ge. (ctrlat + (0.5*dlat))) cycle
            found = .true.
            exit
         end do ! r
         if (.not. found) then
            this%qc(j) = QC_REJECT
            write(LIS_logunit,*) &
                 '[INFO] find_LIS_cols_rows rejecting observation j: ',j,&
                 ' net: ',trim(this%net(j)), &
                 ' platform: ',trim(this%platform(j)), &
                 ' lat: ',this%lat(j), &
                 ' lon: ',this%lon(j)
            cycle
         end if

         ! Then longitude
         found = .false.
         do c = 1, LIS_rc%gnc(nest)
            gindex = c + (r-1)*LIS_rc%gnc(nest)
            ctrlon = LIS_domain(nest)%glon(gindex)

            ! EMK 8 June 2017...Imprecise floating point arithmetic
            ! can cause an ob with a longitude exactly lying on a grid box
            ! edge to "fall between the cracks".  To avoid this, we generally
            ! only check if the ob is west of the eastern grid box edge.
            ! An exception is along the first column, where we need to see
            ! if the ob is west of the LIS domain.
            if (c .eq. 1) then
               if (this%lon(j) .lt. (ctrlon - (0.5*dlon))) cycle
            end if
            if (this%lon(j) .ge. (ctrlon + (0.5*dlon))) cycle
            found = .true.
            exit
         end do ! c
         if (.not. found) then
            this%qc(j) = QC_REJECT
            write(LIS_logunit,*) &
                 '[INFO] find_LIS_cols_rows rejecting observation j: ',j,&
                 ' net: ',trim(this%net(j)), &
                 ' platform: ',trim(this%platform(j)), &
                 ' lat: ',this%lat(j), &
                 ' lon: ',this%lon(j)
            cycle
         end if

         ! Save the column/row for sharing with other PETs
         cols_pet(j) = c
         rows_pet(j) = r

      end do ! j

      ! Collect cols
      allocate(cols(nobs))
      cols = 0
#if (defined SPMD)
      call MPI_Barrier(LIS_MPI_COMM, ierr)
      call handle_mpi_error(ierr, &
           'MPI_Barrier call in find_LIS_cols_rows')
      call MPI_ALLREDUCE(cols_pet,cols,nobs,MPI_INTEGER, &
           MPI_SUM, LIS_MPI_COMM, ierr)
      call handle_mpi_error(ierr, &
           'MPI_Allreduce call in find_LIS_cols_rows')
#else
      cols(:) = cols_pet(:)
#endif
      deallocate(cols_pet)

      ! Collect rows
      allocate(rows(nobs))
      rows = 0
#if (defined SPMD)
      call MPI_Barrier(LIS_MPI_COMM, ierr)
      call handle_mpi_error(ierr, &
           'MPI_Barrier call in find_LIS_cols_rows')
      call MPI_ALLREDUCE(rows_pet,rows,nobs,MPI_INTEGER, &
           MPI_SUM, LIS_MPI_COMM, ierr)
      call handle_mpi_error(ierr, &
           'MPI_Allreduce call in find_LIS_cols_rows')
#else
      rows(:) = rows_pet(:)
#endif
      deallocate(rows_pet)

   end subroutine find_LIS_cols_rows

   !---------------------------------------------------------------------------
   ! Reject observations over water
   subroutine USAF_waterQC(this,nest,silent_rejects)

      ! Imports
      use LIS_LMLCMod, only: LIS_LMLC
      use LIS_logMod, only: LIS_logunit
      use LIS_mpiMod

      ! Defaults
      implicit none

      ! Arguments
      type(USAF_ObsData), intent(inout) :: this
      integer, intent(in) :: nest
      logical,optional,intent(in) :: silent_rejects

      ! Local variables
      integer, allocatable :: cols(:)
      integer, allocatable :: rows(:)
      integer :: nobs
      integer :: j
      integer :: ierr
      double precision :: t1, t2
      integer :: reject_count
      logical :: silent_rejects_local

      ! Sanity check
      nobs = this%nobs
      if (nobs .eq. 0) then
         write(LIS_logunit,*)&
              '[INFO] waterQC found no observations to test'
         return
      end if

      silent_rejects_local = .false.
      if (present(silent_rejects)) then
         silent_rejects_local = silent_rejects
      end if

#if (defined SPMD)
      call MPI_Barrier(LIS_MPI_COMM, ierr)
      call handle_mpi_error(ierr, &
           'MPI_Barrier call in USAF_waterQC')
      t1 = MPI_Wtime()
#endif

      ! Find the LIS column and row for each observation
      call find_LIS_cols_rows(this,nest,cols,rows)

      ! Screen out the observations that are over water.
      ! TODO: Parallelize.
      reject_count = 0
      do j = 1, nobs

         ! Observation not in LIS domain
         if (cols(j) .eq. 0) cycle
         if (rows(j) .eq. 0) cycle

         ! Skip if this is a land point
         if (LIS_LMLC(nest)%glandmask(cols(j),rows(j)) .gt. 0) cycle

         this%qc(j) = QC_REJECT
         reject_count = reject_count + 1

         if (.not. silent_rejects_local) then
            write(LIS_logunit,*) &
                 '[INFO] waterQC rejecting observation i: ',j, &
                 ' net: ',trim(this%net(j)), &
                 ' platform: ',trim(this%platform(j)), &
                 ' lat: ',this%lat(j), &
                 ' lon: ',this%lon(j), &
                 ' obs: ',this%obs(j)
         end if
      end do ! j

      write(LIS_logunit,*)&
              '[INFO] waterQC rejected ',reject_count,' observations'

#if (defined SPMD)
      call MPI_Barrier(LIS_MPI_COMM, ierr)
      call handle_mpi_error(ierr, &
           'MPI_Barrier call in USAF_waterQC')
      t2 = MPI_Wtime()
      write(LIS_logunit,*) &
           '[INFO] Elapsed time in waterQC is ',t2-t1,' seconds'
#endif

      ! Clean up
      if (allocated(cols)) deallocate(cols)
      if (allocated(rows)) deallocate(rows)

   end subroutine USAF_waterQC

   !---------------------------------------------------------------------------
   ! Temperature check to identify potential snowfall.  Based on Lopez (2013)
   ! and earlier operational AGRMET algorithm.
   subroutine USAF_snowQC(this,nest,hourindex,threshold,silent_rejects)

      ! Imports
      use AGRMET_forcingMod, only: agrmet_struc
      use LIS_coreMod, only: LIS_rc, LIS_ews_halo_ind, LIS_ewe_halo_ind, &
           LIS_nss_halo_ind, LIS_nse_halo_ind, LIS_localPet, &
           LIS_ews_ind, LIS_ewe_ind, LIS_nss_ind, LIS_nse_ind
      use LIS_logMod, only: LIS_logunit
      use LIS_mpiMod

      ! Defaults
      implicit none

      ! Arguments
      type(USAF_ObsData), intent(inout) :: this
      integer, intent(in) :: nest
      integer, intent(in) :: hourindex
      real,optional,intent(in) :: threshold
      logical,optional,intent(in) :: silent_rejects

      ! Local variables
      type(hash_list), allocatable :: hash2d(:,:)
      real, allocatable :: sfctmp_pet(:,:), sfctmp(:,:)
      real, allocatable :: sfctmp_1d_pet(:), sfctmp_1d(:)
      integer :: cmax,rmax
      integer :: nobs
      integer :: r,c,j
      integer :: rstart,rend,cstart,cend
      integer :: ierr
      double precision :: t1, t2
      integer :: nobs_cr,job
      integer, allocatable :: jobs_cr_vector(:)
      integer :: reject_count
      real :: threshold_local
      logical :: silent_rejects_local

      ! Sanity check
      nobs = this%nobs
      if (nobs .eq. 0) then
         write(LIS_logunit,*)&
              '[INFO] snowQC found no observations to test'
         return
      end if

#if (defined SPMD)
      call MPI_Barrier(LIS_MPI_COMM, ierr)
      call handle_mpi_error(ierr, &
           'MPI_Barrier call in USAF_snowQC')
      t1 = MPI_Wtime()
#endif

      ! Here we create a 2d hash table storing the index values of each ob
      ! in linked lists for each LIS grid box.
      call build_hash2d(this,nest,cmax,rmax,hash2d)

      ! Find the local (non-halo) bounds in the global grid
      cstart = 1 - &
           (LIS_ews_halo_ind(nest,LIS_localPet+1) - &
            LIS_ews_ind(nest,LIS_localPet+1))
      cend = LIS_rc%lnc(nest) - &
           (LIS_ewe_halo_ind(nest,LIS_localPet+1) - &
            LIS_ewe_ind(nest,LIS_localPet+1))
      rstart = 1 - &
           (LIS_nss_halo_ind(nest,LIS_localPet+1) - &
            LIS_nss_ind(nest,LIS_localPet+1))
      rend = LIS_rc%lnr(nest) - &
           (LIS_nse_halo_ind(nest,LIS_localPet+1) - &
            LIS_nse_ind(nest,LIS_localPet+1))

      ! Collect the global surface temperature analysis
      allocate(sfctmp_pet(LIS_rc%gnc(nest), LIS_rc%gnr(nest)))
      sfctmp_pet(:,:) = 0
      sfctmp_pet( &
           LIS_ews_ind(nest,LIS_localPet+1):LIS_ewe_ind(nest,LIS_localPet+1), &
           LIS_nss_ind(nest,LIS_localPet+1):LIS_nse_ind(nest,LIS_localPet+1)) &
           = &
           agrmet_struc(nest)%sfctmp(hourindex,cstart:cend,rstart:rend)
      allocate(sfctmp_1d_pet(LIS_rc%gnc(nest)*LIS_rc%gnr(nest)))
      sfctmp_1d_pet(:) = 0
      do r = 1, LIS_rc%gnr(nest)
         do c = 1, LIS_rc%gnc(nest)
            sfctmp_1d_pet(c + (r-1)*LIS_rc%gnc(nest)) = &
                 sfctmp_pet(c,r)
         end do ! c
      end do ! r
      deallocate(sfctmp_pet)

#if (defined SPMD)
      allocate(sfctmp_1d(LIS_rc%gnc(nest)*LIS_rc%gnr(nest)))
      sfctmp_1d(:) = 0
      call MPI_Barrier(LIS_MPI_COMM, ierr)
      call handle_mpi_error(ierr, &
           'MPI_Barrier call in USAF_snowQC')
      call MPI_ALLREDUCE(sfctmp_1d_pet,sfctmp_1d, &
           LIS_rc%gnc(nest)*LIS_rc%gnr(nest),MPI_REAL, &
           MPI_SUM, LIS_MPI_COMM, ierr)
      call handle_mpi_error(ierr, &
           'MPI_Allreduce call in USAF_snowQC')
      deallocate(sfctmp_1d_pet)
#endif
      allocate(sfctmp(LIS_rc%gnc(nest), LIS_rc%gnr(nest)))
      sfctmp(:,:) = 0
      do r = 1, LIS_rc%gnr(nest)
         do c = 1, LIS_rc%gnc(nest)
            sfctmp(c,r) = &
                 sfctmp_1d(c + (r-1)*LIS_rc%gnc(nest))
         end do ! c
      end do ! r
      deallocate(sfctmp_1d)

      ! Now set the temperature threshold
      threshold_local = 275.
      if (present(threshold)) then
         threshold_local = threshold
      end if

      silent_rejects_local = .false.
      if (present(silent_rejects)) then
         silent_rejects_local = silent_rejects
      end if

      ! Now loop through the global grid, identifying snow points
      reject_count = 0
      do r = 1, LIS_rc%gnr(nest)
         do c = 1, LIS_rc%gnc(nest)

            ! Make sure obs are actually in this box.
            if (hash2d(c,r)%obindex .eq. MISSING) cycle

            ! Skip if land mask is water -- in practice, no analyzed surface
            ! temperature is available.
!            if (LIS_LMLC(nest)%glandmask(c,r) .le. 0) cycle

            ! See if this is a snowy location
            if (sfctmp(c,r) .gt. threshold_local) cycle

            ! Get list of obs in current grid box
            call get_hash2d(LIS_rc%gnc(nest),LIS_rc%gnr(nest),hash2d,c,r,&
                 nobs_cr,jobs_cr_vector)
            if (nobs_cr .eq. 0) cycle

            do j = 1, nobs_cr
               job = jobs_cr_vector(j)
               this%qc(job) = QC_REJECT
               reject_count = reject_count + 1

               if (.not. silent_rejects_local) then
                  write(LIS_logunit,*) &
                       '[INFO] snowQC rejecting observation i: ',job, &
                       ' net: ',trim(this%net(job)), &
                       ' platform: ',trim(this%platform(job)), &
                       ' lat: ',this%lat(job), &
                       ' lon: ',this%lon(job), &
                       ' obs: ',this%obs(job), &
                       ' sfcT: ',sfctmp(c,r)
               end if
            end do ! j

            deallocate(jobs_cr_vector)
         end do ! c
      end do ! r

      ! Clean up
      call destroy_hash2d(LIS_rc%gnc(nest),LIS_rc%gnr(nest),hash2d)
      deallocate(hash2d)
      deallocate(sfctmp)

      if (reject_count .gt. 0) then
         write(LIS_logunit,*)&
              '[INFO] snowQC rejected ',reject_count,' observations'
      end if

#if (defined SPMD)
      call MPI_Barrier(LIS_MPI_COMM, ierr)
      call handle_mpi_error(ierr, &
           'MPI_Barrier call in USAF_snowQC')
      t2 = MPI_Wtime()
      write(LIS_logunit,*) &
           '[INFO] Elapsed time in snowQC is ',t2-t1,' seconds'
#endif

   end subroutine USAF_snowQC

   !---------------------------------------------------------------------------
   ! Checks for snow on ground.  Based on earlier operational AGRMET code.
   !
   ! This implementation is parallelized by distributing work among the
   ! different LIS ESMF PETs.  To improve performance, a 2D hash table is
   ! constructed to group observations by LIS grid box; this allows for
   ! faster comparison with LIS snow depth.
   subroutine USAF_snowDepthQC(this,nest,silent_rejects)

      ! Imports
      use LIS_coreMod, only: LIS_domain, LIS_rc, &
           LIS_ews_halo_ind, LIS_ewe_halo_ind, &
           LIS_nss_halo_ind, LIS_nse_halo_ind, &
           LIS_localPet, &
           LIS_ews_ind, LIS_ewe_ind, LIS_nss_ind, LIS_nse_ind
      use LIS_logMod, only: LIS_logunit
      use LIS_snowMod, only: LIS_snow_struc
      use LIS_mpiMod

      ! Defaults
      implicit none

      ! Arguments
      type(USAF_ObsData), intent(inout) :: this
      integer, intent(in) :: nest
      logical,optional,intent(in) :: silent_rejects

      ! Local variables
      type(hash_list), allocatable :: hash2d(:,:)
      real, allocatable :: snowdepth_pet(:,:), snowdepth(:,:)
      real, allocatable :: snowdepth_1d_pet(:), snowdepth_1d(:)
      integer :: cmax,rmax
      integer :: nobs
      integer :: r,c,j
      integer :: rstart,rend,cstart,cend
      integer :: ierr
      double precision :: t1, t2
      integer :: nobs_cr,job
      integer, allocatable :: jobs_cr_vector(:)
      integer :: reject_count
      logical :: silent_rejects_local
      integer :: rglb,cglb
      integer :: gid

      ! Sanity check
      nobs = this%nobs
      if (nobs .eq. 0) then
         write(LIS_logunit,*)&
              '[INFO] snowDepthQC found no observations to test'
         return
      endif

#if (defined SPMD)
      call MPI_Barrier(LIS_MPI_COMM, ierr)
      call handle_mpi_error(ierr, &
           'MPI_Barrier call in USAF_snowDepthQC')
      t1 = MPI_Wtime()
#endif

      ! Here we create a 2d hash table storing the index values of each ob
      ! in linked lists for each LIS grid box.
      call build_hash2d(this,nest,cmax,rmax,hash2d)

      ! Find the local (non-halo) bounds in the global grid
      cstart = 1 - &
           (LIS_ews_halo_ind(nest,LIS_localPet+1) - &
            LIS_ews_ind(nest,LIS_localPet+1))
      cend = LIS_rc%lnc(nest) - &
           (LIS_ewe_halo_ind(nest,LIS_localPet+1) - &
            LIS_ewe_ind(nest,LIS_localPet+1))
      rstart = 1 - &
           (LIS_nss_halo_ind(nest,LIS_localPet+1) - &
            LIS_nss_ind(nest,LIS_localPet+1))
      rend = LIS_rc%lnr(nest) - &
           (LIS_nse_halo_ind(nest,LIS_localPet+1) - &
            LIS_nse_ind(nest,LIS_localPet+1))

      ! Collect the global snow depth field
      allocate(snowdepth_pet(LIS_rc%gnc(nest), LIS_rc%gnr(nest)))
      snowdepth_pet(:,:) = 0
      rglb = LIS_nss_ind(nest,LIS_localPet+1)
      do r = rstart,rend
         cglb = LIS_ews_ind(nest,LIS_localPet+1)
         do c = cstart,cend
            gid = LIS_domain(nest)%gindex(c,r)
            if (gid .ne. -1) then
               snowdepth_pet(cglb,rglb) = LIS_snow_struc(nest)%snowdepth(gid)
            end if
            cglb = cglb + 1
         end do ! c
         rglb = rglb + 1
      end do ! r

      allocate(snowdepth_1d_pet(LIS_rc%gnc(nest)*LIS_rc%gnr(nest)))
      snowdepth_1d_pet(:) = 0
      do r = 1, LIS_rc%gnr(nest)
         do c = 1, LIS_rc%gnc(nest)
            snowdepth_1d_pet(c + (r-1)*LIS_rc%gnc(nest)) = &
                 snowdepth_pet(c,r)
         end do ! c
      end do ! r
      deallocate(snowdepth_pet)

#if (defined SPMD)
      allocate(snowdepth_1d(LIS_rc%gnc(nest)*LIS_rc%gnr(nest)))
      snowdepth_1d(:) = 0
      call MPI_Barrier(LIS_MPI_COMM, ierr)
      call handle_mpi_error(ierr, &
           'MPI_Barrier call in USAF_snowDepthQC')
      call MPI_ALLREDUCE(snowdepth_1d_pet,snowdepth_1d, &
           LIS_rc%gnc(nest)*LIS_rc%gnr(nest),MPI_REAL, &
           MPI_SUM, LIS_MPI_COMM, ierr)
      call handle_mpi_error(ierr, &
           'MPI_Allreduce call in USAF_snowDepthQC')
      deallocate(snowdepth_1d_pet)
#endif
      allocate(snowdepth(LIS_rc%gnc(nest), LIS_rc%gnr(nest)))
      snowdepth(:,:) = 0
      do r = 1, LIS_rc%gnr(nest)
         do c = 1, LIS_rc%gnc(nest)
            snowdepth(c,r) = &
                 snowdepth_1d(c + (r-1)*LIS_rc%gnc(nest))
         end do ! c
      end do ! r
      deallocate(snowdepth_1d)

      silent_rejects_local = .false.
      if (present(silent_rejects)) then
         silent_rejects_local = silent_rejects
      end if

      ! Now loop through the global grid, identifying snow points
      reject_count = 0
      do r = 1, LIS_rc%gnr(nest)
         do c = 1, LIS_rc%gnc(nest)

            ! Make sure obs are actually in this box.
            if (hash2d(c,r)%obindex .eq. MISSING) cycle

            ! Skip if land mask is water
!            if (LIS_LMLC(nest)%glandmask(c,r) .le. 0) cycle

            ! See if this is a snowy location
            if (snowdepth(c,r) .le. 0) cycle

            ! Get list of obs in current grid box
            call get_hash2d(LIS_rc%gnc(nest),LIS_rc%gnr(nest),hash2d,c,r,&
                 nobs_cr,jobs_cr_vector)
            if (nobs_cr .eq. 0) cycle

            do j = 1, nobs_cr
               job = jobs_cr_vector(j)
               this%qc(job) = QC_REJECT
               reject_count = reject_count + 1

               if (.not. silent_rejects_local) then
                  write(LIS_logunit,*) &
                       '[INFO] snowDepthQC rejecting observation i: ',job, &
                       ' net: ',trim(this%net(job)), &
                       ' platform: ',trim(this%platform(job)), &
                       ' lat: ',this%lat(job), &
                       ' lon: ',this%lon(job), &
                       ' obs: ',this%obs(job), &
                       ' snowdepth: ',snowdepth(c,r)
               end if
            end do ! j

            deallocate(jobs_cr_vector)
         end do ! c
      end do ! r

      ! Clean up
      call destroy_hash2d(LIS_rc%gnc(nest),LIS_rc%gnr(nest),hash2d)
      deallocate(hash2d)
      deallocate(snowdepth)

      write(LIS_logunit,*)'[INFO] snowDepthQC rejected ',reject_count,&
           ' observations'

#if (defined SPMD)
      call MPI_Barrier(LIS_MPI_COMM, ierr)
      call handle_mpi_error(ierr, &
           'MPI_Barrier call in USAF_snowDepthQC')
      t2 = MPI_Wtime()
      write(LIS_logunit,*) &
           '[INFO] Elapsed time in snowDepthQC is ',t2-t1,' seconds'
#endif

   end subroutine USAF_snowDepthQC

   !---------------------------------------------------------------------------
   ! Reset small non-zero values to zero.
   subroutine zeroTrace(nr,nc,mrgp)
      implicit none
      integer,intent(in) :: nr
      integer,intent(in) :: nc
      real, intent(inout) :: mrgp(nc,nr)
      integer :: r,c
      do r = 1,nr
         do c = 1,nc
            if (mrgp(c,r) .lt. 0.001) then
               mrgp(c,r) = 0
            end if
         end do ! c
      end do ! r
   end subroutine zeroTrace

   !---------------------------------------------------------------------------
   ! Performs analysis on a screen-level variable (e.g., 2-m temperature,
   ! 2-m RH, 10-m wind speed) using a NWP background field and irregularly
   ! positioned observations.
   !
   ! (1) Observations are subjected to quality control tests:
   !     (a) Reports are rejected over water. [waterQC]
   !     (b) Reports are checked for duplicates. [dupQC]
   !     (c) Observations are compared with the background field and are
   !         rejected if the absolute differences are too high, suggesting
   !         gross error. [backQC]
   !     (d) Superobservations are created for observations of the same type
   !         in the same LIS grid box. [superstatQC]
   ! (2) All observations (including superobservations) that passed quality
   !     control are merged into a new ObsData structure.
   ! (3) The inverse data density is calculated for each observation (see
   !     Bratseth 1986).
   ! (4) The Bratseth analysis is generated via iteration at each observation
   !     point.  This successive correction scheme converges to Optimal
   !     Interpolation without direct matrix inversion, providing significant
   !     computational savings and allowing a psuedo-global analysis (the
   !     radius of influence is specified by the background error covariance
   !     rather than arbitrary limits on the number of observations). The
   !     observed, background, and analysis values are also collected in a
   !     OBA structure for later output to file.
   ! (5) The Bratseth analysis is interpolated to the LIS grid in a single
   !     pass, using an algorithm similar to that of Daley (1991),
   !     Pedder (1993), and Kalnay (2003).
   !
   ! NOTES:
   ! (1) This implementation uses a Gaussian function to model error
   !     covariances.  Semivariogram analyses suggest the inverse exponential
   !     function has a better fit to actual statistics, but the Gaussian
   !     function has a much shorter radius of influence that greatly speeds
   !     up the analysis. 
   !---------------------------------------------------------------------------
   subroutine USAF_analyzeScreen(screenObs,nest,back,sigmaBSqr, &
        max_dist, backErrScaleLength,analysis, screenOBA)

      ! Imports
      use AGRMET_forcingMod, only: agrmet_struc
      use LIS_coreMod, only: LIS_rc, LIS_ews_halo_ind, LIS_ewe_halo_ind, &
           LIS_nss_halo_ind, LIS_nse_halo_ind, LIS_localPet
      use LIS_logMod, only: LIS_logunit
      use USAF_OBAMod, only: OBA, createOBA, assignOBA

      ! Defaults
      implicit none

      ! Arguments
      type(USAF_ObsData),intent(inout) :: screenObs
      integer,intent(in) :: nest
      real, intent(in) :: back(LIS_rc%gnc(nest), LIS_rc%gnr(nest))
      real, intent(in) :: sigmaBSqr
      real, intent(in) :: max_dist
      real, intent(in) :: backErrScaleLength
      real, intent(inout) :: analysis(LIS_rc%lnc(nest),LIS_rc%lnr(nest))
      type(OBA), intent(out) :: screenOBA

      ! Local variables
      type(USAF_ObsData) :: screenObsGood
      integer :: goodObs
      integer :: nobs
      real, allocatable :: invDataDensities(:)
      real, allocatable :: sumObsEstimates(:)
      integer :: npasses
      character(len=32) :: new_name,type
      real :: convergeThresh
      integer :: j

      TRACE_ENTER("bratseth_analyzeScrn")
      ! Initialize analysis field with the background first guess.  This will
      ! be changed below as needed.
      analysis(:,:) = back(LIS_ews_halo_ind(nest,LIS_localPet+1): &
           LIS_ewe_halo_ind(nest,LIS_localPet+1), &
           LIS_nss_halo_ind(nest,LIS_localPet+1): &
           LIS_nse_halo_ind(nest,LIS_localPet+1))

      ! If no observations, skip the Bratseth analysis
      if ( screenObs%nobs .eq. 0) then
         if (agrmet_struc(nest)%oba_switch .eq. 1 .or. &
              agrmet_struc(nest)%oba_switch .eq. 2) then
            screenOBA = createOBA(nest,maxobs=0)
         end if
         call USAF_destroyObsData(screenObs)
         TRACE_EXIT("bratseth_analyzeScrn")
         return
      end if

      ! Reject observations over water
      write(LIS_logunit,*)'[INFO] Running waterQC on surface observations'
      call USAF_waterQC(screenObs,nest)

      ! Handle duplicate reports
      write(LIS_logunit,*)'[INFO] Running dupQC on surface observations'
      call USAF_dupQC(screenObs)

      ! EMK...Option 2 just captures O and B info, and skips the
      ! analysis.
      if (agrmet_struc(nest)%oba_switch .eq. 2) then
         goodObs = USAF_countGoodObs(screenObs)
         screenOBA = createOBA(nest,maxobs=goodObs)
         do j = 1, screenObs%nobs
            if (screenObs%qc(j) .eq. QC_REJECT) cycle
            call assignOBA(screenOBA, &
                 screenObs%net(j), screenObs%platform(j), &
                 screenObs%lat(j), screenObs%lon(j), &
                 screenObs%obs(j), screenObs%back(j), &
                 A=0.)
         end do ! j
         call USAF_destroyObsData(screenObs)
         TRACE_EXIT("bratseth_analyzeScrn")
         return ! EMK TEST for O-B
      end if

      ! Compare with background field
      write(LIS_logunit,*)'[INFO] Running backQC on surface observations'
      call USAF_backQC(screenObs,sigmaBSqr)

      ! Create "superobservations" from close reports
      ! FIXME:  Pass network argument.
      new_name = "SUPEROB"
      write(LIS_logunit,*)'[INFO] Running superstatQC on surface observations'
      call USAF_superstatQC(screenObs,nest,new_name)
      type = new_name
      call USAF_interpBackToTypeObsData(screenObs,nest, &
           LIS_rc%gnc(nest),LIS_rc%gnr(nest),back,type)

      ! At this point, QC is done.  Copy the good obs into a new structure
      ! for the analysis (this will speed up analysis calculations by
      ! filtering bad obs from the arrays.)
      goodObs = USAF_countGoodObs(screenObs)
      call USAF_createObsData(screenObsGood,nest,maxobs=goodObs)
      call USAF_filterObsData(screenObsGood,screenObs)
      call USAF_destroyObsData(screenObs)
      nobs = screenObsGood%nobs
      if (nobs .eq. 0) then
         if (agrmet_struc(nest)%oba_switch .eq. 1) then
            screenOBA = createOBA(nest,maxobs=0)
         end if
         call USAF_destroyObsData(screenObsGood)
         TRACE_EXIT("bratseth_analyzeScrn")
         return
      end if

      ! Calculate (inverse) data density around each observation.
      call calc_invDataDensities(screenObsGood,sigmaBSqr,nest, &
           max_dist,backErrScaleLength,is_stn,invDataDensities)

      ! Run Bratseth analysis at observation points, and collect the sum of
      ! the corrections at each observation point (in sumObsEstimates), along
      ! with the required number of iterations (npasses).  Also return
      ! OBA information for output.
      convergeThresh = 0.01
      call calc_obsAnalysis(screenObsGood,sigmaBSqr,nobs,invDataDensities, &
           nest, max_dist, backErrScaleLength,convergeThresh,is_stn,&
           sumObsEstimates, npasses, screenOBA)

      ! Calculate analysis at grid points.
      call calc_gridAnalysis(screenObsGood,nest,sigmaBSqr,nobs, &
           invDataDensities,sumObsEstimates,npasses,back,max_dist,&
           backErrScaleLength,analysis)

      ! Clean up
      deallocate(invDataDensities)
      deallocate(sumObsEstimates)
      call USAF_destroyObsData(screenObsGood)
      TRACE_EXIT("bratseth_analyzeScrn")

   end subroutine USAF_analyzeScreen

   !---------------------------------------------------------------------------
   ! Copy observations to new ObsData structure, filtering out those that have
   ! been flagged for rejection by quality control.
   subroutine USAF_filterObsData(this,obsData)

      ! Defaults
      implicit none

      ! Arguments
      type(USAF_ObsData), intent(inout) :: this
      type(USAF_ObsData), intent(in) :: obsData

      ! Local variables
      integer :: j

      if ( obsData%nobs .eq. 0) return

      do j = 1, obsData%nobs
         if (obsData%qc(j) .eq. QC_REJECT) cycle
         call USAF_assignObsData(this, &
              obsData%net(j),obsData%platform(j), &
              obsData%obs(j),obsData%lat(j),obsData%lon(j), &
              obsData%sigmaOSqr(j), &
              obsData%oErrScaleLength(j), &
              back = obsData%back(j))
      end do ! j

   end subroutine USAF_filterObsData

   !---------------------------------------------------------------------------
   ! Read 8-km CMORPH data and copy to ObsData structures.
   subroutine USAF_getCMORPHObsData(nest,j6hr,use_twelve, &
        precip3, precip6, precip9, precip12, pcp_src)

      ! Imports
      use AGRMET_forcingMod, only : agrmet_struc
      use LIS_logMod, only : LIS_logunit, LIS_getNextUnitNumber, &
           LIS_releaseUnitNumber
      use LIS_timeMgrMod, only : LIS_julhr_date

      ! Defaults
      implicit none

      ! Arguments
      integer, intent(in) :: nest
      integer, intent(in) :: j6hr
      logical, intent(in) :: use_twelve
      type(USAF_obsData), intent(inout) :: precip3
      type(USAF_obsData), intent(inout) :: precip6
      type(USAF_obsData), intent(inout) :: precip9
      type(USAF_obsData), intent(inout) :: precip12
      character(len=6), intent(in) :: pcp_src(4)

      ! Local variables
      integer, parameter :: XD = 4948
      integer, parameter :: YD = 1649
      integer, parameter :: NCMOR = XD*YD
      real :: precip(XD,YD)
      character(len=32) :: net, platform
      real :: sigmaOSqr, oErrScaleLength
      integer :: count_good_obs
      character(len=120) :: fname
      real :: ob
      real :: rlat, rlon
      integer :: yr, mo, da, hr
      integer :: j3hr
      integer :: ftn, ios
      integer :: i,j,k

      external :: cmorfile_agrmet

      TRACE_ENTER("bratseth_getCMORPH")
      net = "CMORPH"
      platform = "CMORPH"
      sigmaOSqr = agrmet_struc(nest)%bratseth_precip_cmorph_sigma_o_sqr
      oErrScaleLength = &
           agrmet_struc(nest)%bratseth_precip_cmorph_err_scale_length

      if (use_twelve) then
         k = 2
      else
         k = 0
      end if

      ! Loop through time levels
      do j3hr = j6hr+3, j6hr+6, 3
         k = k + 1

         ! Set Bratseth error statistics based on source of background field.
         call USAF_setBratsethPrecipStats(pcp_src(k),nest)

         if (agrmet_struc(nest)%cmorswch .eq. 1) then

            ! Get CMORPH file name
            call LIS_julhr_date( j3hr, yr,mo,da,hr)
            call cmorfile_agrmet( fname, agrmet_struc(nest)%agrmetdir, &
                 agrmet_struc(nest)%cmordir, &
                 agrmet_struc(nest)%use_timestamp, &
                 yr, mo, da, hr )

            precip(:,:) = MISSING
            count_good_obs = 0

            ! Read the data
            ftn = LIS_getNextUnitNumber()
            open(unit=ftn,file=fname, status='old',access='direct', &
                 form='unformatted',recl=XD*YD*4,iostat=ios)
            if (ios .eq. 0) then
               read (ftn,rec=1) precip
               ! Note:  i is latitude, j is longitude, file is written
               ! longitude, latitude.
               ! Note:  Northern ten rows do not have trustworthy data based
               ! on examination of JJA 2012 files.
!               do j = 1,YD
               do j = 1,YD-10
                  do i = 1,XD
                     if (precip(i,j) .lt. 0) then
                        precip(i,j)=MISSING
                        cycle
                     end if

                     ! Get the lat and lon of good CMORPH ob.
                     ! Logic based on CPC GrADS control file for 8km CMORPH

                     ! EMK 20120809...Examination of output suggests the
                     ! CMORPH data are shifted south and west by about 10
                     ! rows, which would explain the bad data along the
                     ! northern edge of the CMORPH domain.  The logic for
                     ! calculating CMORPH lat and long is therefore adjusted
                     ! to compensate.

                     ! Original code
!                     rlat = -59.963614 + (j-1)*0.072771377
!                     rlon = 0.036378335 + (i-1)*0.072756669
!                     if (rlon > 180.) rlon = rlon - 360.

                     ! To fix phase shift
                     rlat = -59.963614 + (j-1+10)*0.072771377
                     rlon = 0.036378335 + (i-1+10)*0.072756669
                     if (rlon > 180.) rlon = rlon - 360.

                     ob = precip(i,j)
                     count_good_obs = count_good_obs + 1

                     ! Now save the ob in the appropriate structure
                     if (k .eq. 1) then
                        call USAF_assignObsData(precip3,net,platform,ob, &
                             rlat,rlon,sigmaOSqr, &
                             oErrScaleLength)
                     else if (k .eq. 2) then
                        call USAF_assignObsData(precip6,net,platform,ob, &
                             rlat,rlon,sigmaOSqr, &
                             oErrScaleLength)
                     else if (k .eq. 3) then
                        call USAF_assignObsData(precip9,net,platform,ob, &
                             rlat,rlon,sigmaOSqr, &
                             oErrScaleLength)
                     else if (k .eq. 4) then
                        call USAF_assignObsData(precip12,net,platform,ob, &
                             rlat,rlon,sigmaOSqr, &
                             oErrScaleLength)
                     end if
                  enddo ! j
               enddo ! i

               write(LIS_logunit,*)'[INFO] Read ',count_good_obs,&
                    ' CMORPH observations'

            else
               write(LIS_logunit,*) &
                    "[ERR] Missing CMORPH precipitation data ", fname
            end if

            call LIS_releaseUnitNumber(ftn)

         end if
      end do ! j3hr
      TRACE_EXIT("bratseth_getCMORPH")

   end subroutine USAF_getCMORPHObsData

   ! Internal subroutine for fetching precipitation from GFS GRIB file.
   ! Based on AGRMET_fldbld_read_precip_gfs.
   ! NOTE:  This assumes the GFS GRIB files have a single precipitation
   ! field that is accumulated from the previous 3 or 6-hours.
   subroutine fldbld_read_precip_gfs( fg_filename, ifguess, jfguess,&
        fg_prec, rc )

      ! Imports
      use LIS_logMod, only : LIS_logunit, LIS_abort, LIS_alert, LIS_verify
#if (defined USE_GRIBAPI)
      use grib_api
#endif

      ! Defaults
      implicit none

      ! Arguments
      character(len=*),  intent(in) :: fg_filename
      integer,        intent(in)    :: ifguess
      integer,        intent(in)    :: jfguess
      real,allocatable,intent(out)   :: fg_prec     (:,:)
      integer, intent(out) :: rc

      ! Locals
      integer                       :: count_prec
      integer                       :: ierr
      integer                       :: k, c, r
      integer                       :: ftn, igrib, nvars
      integer                       :: pds7_val, pds8_val, pds9_val
      ! EMK...For GRIB 2
      integer :: editionNumber_val
      integer :: param_disc_val, param_cat_val, &
           param_num_val
      real,           allocatable   :: dum1d       ( : )
      logical :: found_inq

      rc = 0
      count_prec = 0

      ! EMK...Before using ECCODES/GRIB_API, see if the GRIB file exists
      ! using a simple inquire statement.  This avoids ECCODES/GRIB_API
      ! writing error messages to stdout/stderr, which may lead to runtime
      ! problems.
      inquire(file=trim(fg_filename),exist=found_inq)
      if (.not. found_inq) then
         write(LIS_logunit,*)'[WARN] Cannot find file '//trim(fg_filename)
         rc = 1
         return
      end if

#if (defined USE_GRIBAPI)
      ! Open the file.
      call grib_open_file(ftn,trim(fg_filename),'r',ierr)
      if (ierr .ne. 0) then
         write(LIS_logunit,*) &
              '[WARN] Failed to open ',trim(fg_filename)
         rc = 1
         return
      end if

      ! Get number of fields in GRIB file.
      call grib_count_in_file(ftn,nvars,ierr)
      if (ierr .ne. 0) then
         write(LIS_logunit,*) &
              '[WARN] error in grib_count_in_file in fldbld_read_precip_gfs'
         call grib_close_file(ftn)
         rc = 1
         return
      end if

      ! Loop through the fields in the file until we find the precipitation.
      ! Error returned from the GRIB code will be interpreted as a corrupt
      ! file, and we will return from the subroutine with rc=1 (after closing
      ! the file).
      do k = 1,nvars

         ! First, get next field record
         call grib_new_from_file(ftn,igrib,ierr)
         if (ierr .ne. 0) then
            write(LIS_logunit,*) &
                 '[WARN] error in grib_new_from_file in fldbld_read_precip_gfs'
            call grib_close_file(ftn)
            rc = 1
            return
         end if

         ! Determine if GRIB1 or GRIB2
         call grib_get(igrib,'editionNumber',editionNumber_val,ierr)
         if (ierr .ne. 0) then
            write(LIS_logunit,*) &
                 '[WARN] error in grib_get in fldbld_read_precip_gfs'
            call grib_release(igrib,ierr)
            call grib_close_file(ftn)
            rc = 1
            return
         end if

         ! Determine if this is the correct field.  Different parameters must
         ! be checked depending on whether this is GRIB1 or GRIB2.
         if (editionNumber_val .eq. 1) then

            call grib_get(igrib,'indicatorOfParameter',pds7_val,ierr)
            if (ierr .ne. 0) then
               write(LIS_logunit,*) &
                    '[WARN] error in grib_get: indicatorOfParameter in ' //&
                    'fldbld_read_precip_gfs'
               call grib_release(igrib,ierr)
               call grib_close_file(ftn)
               rc = 1
               return
            end if
            if (pds7_val .ne. 61) then
               call grib_release(igrib,ierr)
               if (ierr .ne. 0) then
                  write(LIS_logunit,*) &
                       '[WARN] error in grib_release in ' //&
                       'fldbld_read_precip_gfs'
                  call grib_close_file(ftn)
                  rc = 1
                  return
               end if
               cycle ! Go to next field in GRIB file
            end if

            call grib_get(igrib,'level',pds9_val,ierr)
            if (ierr .ne. 0) then
               write(LIS_logunit,*) &
                    '[WARN] error in grib_get: level in ' //&
                    'fldbld_read_precip_gfs'
               call grib_release(igrib,ierr)
               call grib_close_file(ftn)
               rc = 1
               return
            end if
            if (pds9_val .ne. 0) then
               call grib_release(igrib,ierr)
               if (ierr .ne. 0) then
                  write(LIS_logunit,*) &
                       '[WARN] error in grib_release in ' //&
                       'fldbld_read_precip_gfs'
                  call grib_close_file(ftn)
                  rc = 1
                  return
               end if
               cycle ! Go to next field in GRIB file
            end if

            call grib_get(igrib,'indicatorOfTypeOfLevel',pds8_val,ierr)
            if (ierr .ne. 0) then
               write(LIS_logunit,*) &
                    '[WARN] error in grib_get: indicatorOfTypeOfLevel in ' //&
                    'fldbld_read_precip_gfs'
               call grib_release(igrib,ierr)
               call grib_close_file(ftn)
               rc = 1
               return
            end if
            if (pds8_val .ne. 1) then
               call grib_release(igrib,ierr)
               if (ierr .ne. 0) then
                  write(LIS_logunit,*) &
                       '[WARN] error in grib_release in ' //&
                       'fldbld_read_precip_gfs'
                  call grib_close_file(ftn)
                  rc = 1
                  return
               end if
               cycle ! Go to next field in GRIB file
            end if

         else ! GRIB2

            call grib_get(igrib,'discipline',param_disc_val,ierr)
            if (ierr .ne. 0) then
               write(LIS_logunit,*) &
                    '[WARN] error in grib_get: discipline in ' //&
                    'fldbld_read_precip_gfs'
               call grib_release(igrib,ierr)
               call grib_close_file(ftn)
               rc = 1
               return
            end if
            if (param_disc_val .ne. 0) then
               call grib_release(igrib,ierr)
               if (ierr .ne. 0) then
                  write(LIS_logunit,*) &
                       '[WARN] error in grib_release in ' //&
                       'fldbld_read_precip_gfs'
                  call grib_close_file(ftn)
                  rc = 1
                  return
               end if
               cycle ! Go to next field in GRIB file
            end if

            call grib_get(igrib,'parameterCategory',param_cat_val,ierr)
            if (ierr .ne. 0) then
               write(LIS_logunit,*) &
                    '[WARN] error in grib_get: parameterCategory in ' //&
                    'fldbld_read_precip_gfs'
               call grib_release(igrib,ierr)
               call grib_close_file(ftn)
               rc = 1
               return
            end if
            if (param_cat_val .ne. 1) then
               call grib_release(igrib,ierr)
               if (ierr .ne. 0) then
                  write(LIS_logunit,*) &
                       '[WARN] error in grib_release in ' //&
                       'fldbld_read_precip_gfs'
                  call grib_close_file(ftn)
                  rc = 1
                  return
               end if
               cycle ! Go to next field in GRIB file
            end if

            call grib_get(igrib,'parameterNumber',param_num_val,ierr)
            if (ierr .ne. 0) then
               write(LIS_logunit,*) &
                    '[WARN] error in grib_get: parameterNumber in ' //&
                    'fldbld_read_precip_gfs'
               call grib_release(igrib,ierr)
               call grib_close_file(ftn)
               rc = 1
               return
            end if
            if (param_num_val .ne. 8) then
               call grib_release(igrib,ierr)
               if (ierr .ne. 0) then
                  write(LIS_logunit,*) &
                       '[WARN] error in grib_release in ' //&
                       'fldbld_read_precip_gfs'
                  call grib_close_file(ftn)
                  rc = 1
                  return
               end if
               cycle ! Go to next field in GRIB file
            end if

         end if ! GRIB vs GRIB2

         ! We found the right field
         allocate ( dum1d   (ifguess*jfguess) )
         call grib_get(igrib, 'values',dum1d,ierr)
         if (ierr .ne. 0) then
            write(LIS_logunit,*) &
                 '[WARN] error in grib_get: values in ' //&
                 'fldbld_read_precip_gfs'
            deallocate(dum1d)
            call grib_release(igrib,ierr)
            call grib_close_file(ftn)
            rc = 1
            return
         end if

         ! At this stage, we have the values of the field.
         call grib_release(igrib,ierr) ! Ignore error here since we have data
         call grib_close_file(ftn)
         allocate(fg_prec(ifguess,jfguess))
         do r=1,jfguess
            do c=1,ifguess
               fg_prec(c,r) = dum1d(c+(r-1)*ifguess)
            enddo
         enddo
         count_prec = count_prec + 1
         deallocate(dum1d)
         exit ! Get out of loop

      end do ! k
#endif

      ! Handle failure to find field
      if (count_prec .eq. 0) then
         rc = 1
         return
      end if

   end subroutine fldbld_read_precip_gfs

   ! Internal subroutine
   subroutine fldbld_read_precip_galwem(fg_filename, ifguess, jfguess,&
        fc_hr, fg_prec, rc )

      ! Imports
      use LIS_logMod, only : LIS_logunit, LIS_abort, LIS_alert, LIS_verify

#if (defined USE_GRIBAPI)
      use grib_api
#endif

      ! Defaults
      implicit none

      ! Arguments
      character(len=*),  intent(in) :: fg_filename
      integer,        intent(in)    :: ifguess
      integer,        intent(in)    :: jfguess
      integer,        intent(in)    :: fc_hr
      real,allocatable,intent(out)   :: fg_prec(:,:)
      integer, intent(out) :: rc

      ! Locals
      integer                       :: count_prec
      integer                       :: ierr
      integer                       :: k
      integer                       :: ftn, igrib, nvars
      integer                       :: param_disc_val, param_cat_val, &
           param_num_val, forecasttime_val
      real,           allocatable   :: dum1d       ( : )
      logical :: found_inq

      rc = 0
      count_prec = 0

      ! EMK...Before using ECCODES/GRIB_API, see if the GRIB file exists
      ! using a simple inquire statement.  This avoids ECCODES/GRIB_API
      ! writing error messages to stdout/stderr, which may lead to runtime
      ! problems.
      inquire(file=trim(fg_filename),exist=found_inq)
      if (.not. found_inq) then
         write(LIS_logunit,*)'[WARN] Cannot find file '//trim(fg_filename)
         rc = 1
         return
      end if

#if (defined USE_GRIBAPI)
      ! Open file
      call grib_open_file(ftn,trim(fg_filename),'r',ierr)
      if (ierr .ne. 0) then
         write(LIS_logunit,*) &
              '[WARN] error from grib_open_file for ' //trim(fg_filename)// &
              ' in fldbld_read_precip_galwem'
         rc = 1
         return
      end if

      ! Get number of fields in this file
      call grib_count_in_file(ftn,nvars,ierr)
      if (ierr .ne. 0) then
         write(LIS_logunit,*) &
              '[WARN] error from grib_count_in_file in ' // &
              'fldbld_read_precip_galwem'
         rc = 1
         return
      end if

      ! Loop through the fields in the file until we find the precipitation.
      ! Make sure we find the accumulation starting from 3-hr prior.
      ! Errors returned from GRIB file will be interpreted as a corrupt file,
      ! and we will return from this subroutine with rc = 1.
      do k = 1,nvars
         ! Find next field in the GRIB file
         call grib_new_from_file(ftn,igrib,ierr)
         if (ierr .ne. 0) then
            write(LIS_logunit,*) &
                 '[WARN] error from grib_new_from_file in ' // &
                 'fldbld_read_precip_galwem'
            call grib_close_file(ftn)
            rc = 1
            return
         end if

         ! See if field is a meteorological product
         call grib_get(igrib,'discipline',param_disc_val,ierr)
         if (ierr .ne. 0) then
            write(LIS_logunit,*) &
                    '[WARN] error from grib_get: discipline in ' // &
                    'fldbld_read_precip_galwem'
            call grib_release(igrib,ierr)
            call grib_close_file(ftn)
            rc = 1
            return
         end if
         if (param_disc_val .ne. 0) then ! wrong field
            call grib_release(igrib,ierr)
            if (ierr .ne. 0) then
               write(LIS_logunit,*) &
                    '[WARN] error from grib_release in ' // &
                    'fldbld_read_precip_galwem'
               call grib_close_file(ftn)
               rc = 1
               return
            end if
            cycle ! Go to next field in GRIB file
         end if

         ! See if field is a moisture category product
         call grib_get(igrib,'parameterCategory',param_cat_val,ierr)
         if (ierr .ne. 0) then
            write(LIS_logunit,*) &
                    '[WARN] error from grib_get: parameterCategory in ' // &
                    'fldbld_read_precip_galwem'
            call grib_release(igrib,ierr)
            call grib_close_file(ftn)
            rc = 1
            return
         end if
         if (param_cat_val .ne. 1) then ! wrong field
            call grib_release(igrib,ierr)
            if (ierr .ne. 0) then
               write(LIS_logunit,*) &
                    '[WARN] error from grib_release in ' // &
                    'fldbld_read_precip_galwem'
               call grib_close_file(ftn)
               rc = 1
               return
            end if
            cycle ! Go to next field in GRIB file
         end if

         ! See if field is total precipitation
         call grib_get(igrib,'parameterNumber',param_num_val,ierr)
         if (ierr .ne. 0) then
            write(LIS_logunit,*) &
                    '[WARN] error from grib_get: parameterNumber in ' // &
                    'fldbld_read_precip_galwem'
            call grib_release(igrib,ierr)
            call grib_close_file(ftn)
            rc = 1
            return
         end if
         if (param_num_val .ne. 8) then ! wrong field
            call grib_release(igrib,ierr)
            if (ierr .ne. 0) then
               write(LIS_logunit,*) &
                    '[WARN] error from grib_release in ' // &
                    'fldbld_read_precip_galwem'
               call grib_close_file(ftn)
               rc = 1
               return
            end if
            cycle ! Go to next field in GRIB file
         end if

         ! Get the forecast time, which here means the time at the start
         ! of the precip accumulation.  This needs to be 3-hr prior to
         ! the current (ending) time fc_hr.
         call grib_get(igrib,'forecastTime',forecasttime_val,ierr)
         if (ierr .ne. 0) then
            write(LIS_logunit,*) &
                    '[WARN] error from grib_get: forecastTime in ' // &
                    'fldbld_read_precip_galwem'
            call grib_release(igrib,ierr)
            call grib_close_file(ftn)
            rc = 1
            return
         end if
         if (forecasttime_val .ne. (fc_hr - 3)) then ! wrong field
            call grib_release(igrib,ierr)
            if (ierr .ne. 0) then
               write(LIS_logunit,*) &
                    '[WARN] error from grib_release in ' // &
                    'fldbld_read_precip_galwem'
               call grib_close_file(ftn)
               rc = 1
               return
            end if
            cycle ! Go to next field in GRIB file
         end if

         ! We found the right field
         allocate ( dum1d   (ifguess*jfguess) )
         call grib_get(igrib,'values',dum1d,ierr)
         if (ierr .ne. 0) then
            write(LIS_logunit,*) &
                 '[WARN] error from grib_get: values in ' // &
                 'fldbld_read_precip_galwem'
            deallocate(dum1d)
            call grib_release(igrib,ierr)
            call grib_close_file(ftn)
            rc = 1
            return
         end if

         ! We have the data.  Wrap this up.
         allocate(fg_prec(ifguess,jfguess))
         fg_prec = reshape(dum1d, (/ifguess,jfguess/))
         deallocate(dum1d)
         count_prec = count_prec + 1
         call grib_release(igrib,ierr) ! We have data, so skip error handling
         call grib_close_file(ftn)
         exit ! Get out of loop

      end do ! k

#endif

      ! Handle failure to find precipitation
      if (count_prec == 0) then
         rc = 1
         return
      end if

      ! At this point, we are done
      rc = 0
      return

   end subroutine fldbld_read_precip_galwem

   ! Public subroutine for selecting appropriate set of Bratseth error
   ! statistics based on source of background field.  This should be called
   ! before each invocation of USAF_analyzePrecip, as the background field
   ! source is dependent on available GRIB files
   subroutine USAF_setBratsethPrecipStats(src,nest)

      ! Imports
      use AGRMET_forcingMod, only : agrmet_struc
      use LIS_logMod, only: LIS_endrun, LIS_logunit

      ! Defaults
      implicit none

      ! Arguments
      character(len=6),intent(in) :: src
      integer, intent(in) :: nest

!      TRACE_ENTER("bratseth_setPrcpStats")
      if (trim(src) == "GALWEM") then
         agrmet_struc(nest)%bratseth_precip_back_err_scale_length = &
              agrmet_struc(nest)%galwem_precip_back_err_scale_length
         agrmet_struc(nest)%bratseth_precip_back_sigma_b_sqr = &
              agrmet_struc(nest)%galwem_precip_back_sigma_b_sqr
         agrmet_struc(nest)%bratseth_precip_gauge_sigma_o_sqr = &
              agrmet_struc(nest)%galwem_precip_gauge_sigma_o_sqr
         agrmet_struc(nest)%bratseth_precip_geoprecip_err_scale_length =&
              agrmet_struc(nest)%galwem_precip_geoprecip_err_scale_length
         agrmet_struc(nest)%bratseth_precip_geoprecip_sigma_o_sqr = &
              agrmet_struc(nest)%galwem_precip_geoprecip_sigma_o_sqr
         agrmet_struc(nest)%bratseth_precip_ssmi_err_scale_length = &
              agrmet_struc(nest)%galwem_precip_ssmi_err_scale_length
         agrmet_struc(nest)%bratseth_precip_ssmi_sigma_o_sqr = &
              agrmet_struc(nest)%galwem_precip_ssmi_sigma_o_sqr
         agrmet_struc(nest)%bratseth_precip_cmorph_err_scale_length = &
              agrmet_struc(nest)%galwem_precip_cmorph_err_scale_length
         agrmet_struc(nest)%bratseth_precip_cmorph_sigma_o_sqr = &
              agrmet_struc(nest)%galwem_precip_cmorph_sigma_o_sqr
         agrmet_struc(nest)%bratseth_precip_imerg_err_scale_length = &
              agrmet_struc(nest)%galwem_precip_imerg_err_scale_length
         agrmet_struc(nest)%bratseth_precip_imerg_sigma_o_sqr = &
              agrmet_struc(nest)%galwem_precip_imerg_sigma_o_sqr
         agrmet_struc(nest)%bratseth_precip_max_dist = &
              agrmet_struc(nest)%galwem_precip_max_dist

      else if (trim(src) == "GFS") then
         agrmet_struc(nest)%bratseth_precip_back_err_scale_length = &
              agrmet_struc(nest)%gfs_precip_back_err_scale_length
         agrmet_struc(nest)%bratseth_precip_back_sigma_b_sqr = &
              agrmet_struc(nest)%gfs_precip_back_sigma_b_sqr
         agrmet_struc(nest)%bratseth_precip_gauge_sigma_o_sqr = &
              agrmet_struc(nest)%gfs_precip_gauge_sigma_o_sqr
         agrmet_struc(nest)%bratseth_precip_geoprecip_err_scale_length =&
              agrmet_struc(nest)%gfs_precip_geoprecip_err_scale_length
         agrmet_struc(nest)%bratseth_precip_geoprecip_sigma_o_sqr = &
              agrmet_struc(nest)%gfs_precip_geoprecip_sigma_o_sqr
         agrmet_struc(nest)%bratseth_precip_ssmi_err_scale_length = &
              agrmet_struc(nest)%gfs_precip_ssmi_err_scale_length
         agrmet_struc(nest)%bratseth_precip_ssmi_sigma_o_sqr = &
              agrmet_struc(nest)%gfs_precip_ssmi_sigma_o_sqr
         agrmet_struc(nest)%bratseth_precip_cmorph_err_scale_length = &
              agrmet_struc(nest)%gfs_precip_cmorph_err_scale_length
         agrmet_struc(nest)%bratseth_precip_cmorph_sigma_o_sqr = &
              agrmet_struc(nest)%gfs_precip_cmorph_sigma_o_sqr
         agrmet_struc(nest)%bratseth_precip_imerg_err_scale_length = &
              agrmet_struc(nest)%gfs_precip_imerg_err_scale_length
         agrmet_struc(nest)%bratseth_precip_imerg_sigma_o_sqr = &
              agrmet_struc(nest)%gfs_precip_imerg_sigma_o_sqr
         agrmet_struc(nest)%bratseth_precip_max_dist = &
              agrmet_struc(nest)%gfs_precip_max_dist

      else
         write(LIS_logunit,*) &
              '[ERR] Unknown source of background precipitation!'
         write(LIS_logunit,*) &
              '[ERR] Source is ',trim(src)
         write(LIS_logunit, *) &
              '[ERR] ABORTING....'
         flush(LIS_logunit)
         call LIS_endrun()
      end if
!      TRACE_EXIT("bratseth_setPrcpStats")

   end subroutine USAF_setBratsethPrecipStats

   ! Public subroutine for selecting appropriate set of Bratseth error
   ! statistics based on source of background field.  This should be called
   ! before each invocation of USAF_analyzeScreen, as the background field
   ! source is dependent on available GRIB files
   subroutine USAF_setBratsethScreenStats(src,n)

      ! Imports
      use AGRMET_forcingMod, only : agrmet_struc
      use LIS_logMod, only: LIS_endrun, LIS_logunit

      ! Defaults
      implicit none

      ! Arguments
      character(len=6),intent(in) :: src
      integer, intent(in) :: n

      TRACE_ENTER("bratseth_setScrnStats")
      if (trim(src) == "GALWEM") then
         agrmet_struc(n)%bratseth_t2m_back_err_scale_length = &
              agrmet_struc(n)%galwem_t2m_back_err_scale_length
         agrmet_struc(n)%bratseth_t2m_back_sigma_b_sqr = &
              agrmet_struc(n)%galwem_t2m_back_sigma_b_sqr

         !write(LIS_logunit,*)'EMK: n, GALWEM T2 err = ', n, &
         !     agrmet_struc(n)%bratseth_t2m_back_sigma_b_sqr
         !flush(LIS_logunit)

         agrmet_struc(n)%bratseth_t2m_stn_sigma_o_sqr = &
              agrmet_struc(n)%galwem_t2m_stn_sigma_o_sqr
         agrmet_struc(n)%bratseth_t2m_max_dist = &
              agrmet_struc(n)%galwem_t2m_max_dist

         agrmet_struc(n)%bratseth_rh2m_back_err_scale_length = &
              agrmet_struc(n)%galwem_rh2m_back_err_scale_length
         agrmet_struc(n)%bratseth_rh2m_back_sigma_b_sqr = &
              agrmet_struc(n)%galwem_rh2m_back_sigma_b_sqr
         agrmet_struc(n)%bratseth_rh2m_stn_sigma_o_sqr = &
              agrmet_struc(n)%galwem_rh2m_stn_sigma_o_sqr
         agrmet_struc(n)%bratseth_rh2m_max_dist = &
              agrmet_struc(n)%galwem_rh2m_max_dist

         agrmet_struc(n)%bratseth_spd10m_back_err_scale_length = &
              agrmet_struc(n)%galwem_spd10m_back_err_scale_length
         agrmet_struc(n)%bratseth_spd10m_back_sigma_b_sqr = &
              agrmet_struc(n)%galwem_spd10m_back_sigma_b_sqr
         agrmet_struc(n)%bratseth_spd10m_stn_sigma_o_sqr = &
              agrmet_struc(n)%galwem_spd10m_stn_sigma_o_sqr
         agrmet_struc(n)%bratseth_spd10m_max_dist = &
              agrmet_struc(n)%galwem_spd10m_max_dist

      else if (trim(src) .eq. "GFS") then
         agrmet_struc(n)%bratseth_t2m_back_err_scale_length = &
              agrmet_struc(n)%gfs_t2m_back_err_scale_length
         agrmet_struc(n)%bratseth_t2m_back_sigma_b_sqr = &
              agrmet_struc(n)%gfs_t2m_back_sigma_b_sqr

         !write(LIS_logunit,*)'EMK: n, GFS T2 err = ', n, &
         !     agrmet_struc(n)%bratseth_t2m_back_sigma_b_sqr
         !flush(LIS_Logunit)

         agrmet_struc(n)%bratseth_t2m_stn_sigma_o_sqr = &
              agrmet_struc(n)%gfs_t2m_stn_sigma_o_sqr
         agrmet_struc(n)%bratseth_t2m_max_dist = &
              agrmet_struc(n)%gfs_t2m_max_dist

         agrmet_struc(n)%bratseth_rh2m_back_err_scale_length = &
              agrmet_struc(n)%gfs_rh2m_back_err_scale_length
         agrmet_struc(n)%bratseth_rh2m_back_sigma_b_sqr = &
              agrmet_struc(n)%gfs_rh2m_back_sigma_b_sqr
         agrmet_struc(n)%bratseth_rh2m_stn_sigma_o_sqr = &
              agrmet_struc(n)%gfs_rh2m_stn_sigma_o_sqr
         agrmet_struc(n)%bratseth_rh2m_max_dist = &
              agrmet_struc(n)%gfs_rh2m_max_dist

         agrmet_struc(n)%bratseth_spd10m_back_err_scale_length = &
              agrmet_struc(n)%gfs_spd10m_back_err_scale_length
         agrmet_struc(n)%bratseth_spd10m_back_sigma_b_sqr = &
              agrmet_struc(n)%gfs_spd10m_back_sigma_b_sqr
         agrmet_struc(n)%bratseth_spd10m_stn_sigma_o_sqr = &
              agrmet_struc(n)%gfs_spd10m_stn_sigma_o_sqr
         agrmet_struc(n)%bratseth_spd10m_max_dist = &
              agrmet_struc(n)%gfs_spd10m_max_dist

      else
         write(LIS_logunit,*) &
              '[ERR] Unknown source of background data!'
         write(LIS_logunit,*) &
              '[ERR] Source is ',trim(src)
         write(LIS_logunit, *) &
              '[ERR] ABORTING....'
         flush(LIS_logunit)
         call LIS_endrun()
      end if
      TRACE_EXIT("bratseth_setScrnStats")

   end subroutine USAF_setBratsethScreenStats

   ! Read bias ratio for background field (S2S version, which uses
   ! data in the LDT parameter file)
   subroutine USAF_pcpBackBiasRatio_s2s(n, yyyymmddhh)

     ! Imports
     use AGRMET_forcingMod, only: agrmet_struc
     use LIS_coreMod, only: LIS_rc
     use LIS_logMod, only: LIS_logunit, LIS_endrun, LIS_verify
#if ( defined USE_NETCDF3 || defined USE_NETCDF4 )
     use netcdf
#endif

     ! Defaults
     implicit none

     ! Arguments
     integer, intent(in) :: n
     character(len=10), intent(in) :: yyyymmddhh

     ! Locals
     integer :: imonth
     logical :: file_exists
     integer :: ncid, ppt_ratio_Id

#if ( defined USE_NETCDF3 || defined USE_NETCDF4 )

     read(yyyymmddhh(5:6),'(i2.2)') imonth

     ! Accumulations ending at 00Z first of month are part of the *previous*
     ! month.  So decrement imonth by one in that situation, and reset to 12
     ! (December) at the year change.
     if (agrmet_struc(n)%pcp_back_bias_ratio_month .ne. 0) then
        if (yyyymmddhh(9:10) .eq. '00' .and. &
             yyyymmddhh(7:8) .eq. '01') then
           imonth = imonth - 1
           if (imonth .eq. 0) then
              imonth = 12
           end if
        end if
     end if
     if (imonth .eq. agrmet_struc(n)%pcp_back_bias_ratio_month .and. &
          agrmet_struc(n)%pcp_back_bias_ratio_month .ne. 0) return

     agrmet_struc(n)%pcp_back_bias_ratio_month = imonth

     inquire(file=LIS_rc%paramfile(n), exist=file_exists)

     if (.not. file_exists) then
        write(LIS_logunit,*) '[ERR] Cannot find ', trim(LIS_rc%paramfile(n))
        call LIS_endrun()
     end if

     write(LIS_logunit,*) '[INFO] Reading precip bias ratio for month ', imonth

     call LIS_verify(nf90_open(path=LIS_rc%paramfile(n),&
               mode=NF90_NOWRITE, ncid=ncid), &
               'Error in nf90_open in USAF_pcpBackBiasRatio_s2s')
     call LIS_verify(nf90_inq_varid(ncid, "PPT_ratio", ppt_ratio_Id), &
          'PPT_ratio field not found in LIS param file')
     call LIS_verify(nf90_get_var(ncid, ppt_ratio_Id, &
          agrmet_struc(n)%pcp_back_bias_ratio, &
               start=(/1,1,imonth/), &
               count=(/LIS_rc%gnc(n),LIS_rc%gnr(n),1/)),&
               'n90_get_var failed for PPT_ratio in LIS param file')
     call LIS_verify(nf90_close(ncid),&
          'nf90_close failed in USAF_pcpBackBiasRatio_s2s')
#endif

   end subroutine USAF_pcpBackBiasRatio_s2s

   ! Read bias ratio for background field (NRT version, which uses
   ! data in standalone netCDF file)
   subroutine USAF_pcpBackBiasRatio_nrt(n, yyyymmddhh)

     ! Imports
     use AGRMET_forcingMod, only: agrmet_struc
     use LIS_coreMod, only: LIS_rc!, LIS_localPet, &
          !LIS_ews_halo_ind, LIS_ewe_halo_ind, &
          !LIS_nss_halo_ind, LIS_nse_halo_ind

     use LIS_logMod, only: LIS_logunit, LIS_endrun, LIS_verify
#if ( defined USE_NETCDF3 || defined USE_NETCDF4 )
     use netcdf
#endif

     ! Defaults
     implicit none

     ! Arguments
     integer, intent(in) :: n
     character(len=10), intent(in) :: yyyymmddhh

     ! Locals
     integer :: imonth
     logical :: file_exists
     integer :: ncid, ppt_ratio_Id

#if ( defined USE_NETCDF3 || defined USE_NETCDF4 )

     read(yyyymmddhh(5:6),'(i2.2)') imonth

     ! Accumulations ending at 00Z first of month are part of the
     ! *previous* month.  So decrement imonth by one in that situation,
     ! and reset to 12 (December) at the year change.
     if (agrmet_struc(n)%pcp_back_bias_ratio_month .ne. 0) then
        if (yyyymmddhh(9:10) .eq. '00' .and. &
             yyyymmddhh(7:8) .eq. '01') then
           imonth = imonth - 1
           if (imonth .eq. 0) then
              imonth = 12
           end if
        end if
     end if
     if (imonth .eq. agrmet_struc(n)%pcp_back_bias_ratio_month .and. &
          agrmet_struc(n)%pcp_back_bias_ratio_month .ne. 0) return

     agrmet_struc(n)%pcp_back_bias_ratio_month = imonth

     ! First, update the GFS ratios
     inquire(file=agrmet_struc(n)%gfs_nrt_bias_ratio_file, &
          exist=file_exists)
     if (.not. file_exists) then
        write(LIS_logunit,*) '[ERR] Cannot find ', &
             trim(agrmet_struc(n)%gfs_nrt_bias_ratio_file)
        call LIS_endrun()
     end if
     call LIS_verify( &
          nf90_open(path=agrmet_struc(n)%gfs_nrt_bias_ratio_file, &
          mode=NF90_NOWRITE, ncid=ncid), &
          'Error in nf90_open in USAF_pcpBackBiasRatio_nrt')
     call LIS_verify(nf90_inq_varid(ncid, "biasRatio", ppt_ratio_Id), &
          'biasRatio field not found in bias file')
     call LIS_verify(nf90_get_var(ncid, ppt_ratio_Id, &
          agrmet_struc(n)%gfs_nrt_bias_ratio, &
          start=(/1,1,imonth/), &
          count=(/LIS_rc%gnc(n),LIS_rc%gnr(n),1/)), &
          'n90_get_var failed for biasRatio in LIS param file')
     write(LIS_logunit,*) &
          '[INFO] Read GFS bias ratios for month ', imonth, &
          ' from ', trim(agrmet_struc(n)%gfs_nrt_bias_ratio_file)
     call LIS_verify(nf90_close(ncid),&
          'nf90_close failed in USAF_pcpBackBiasRatio_nrt')

     ! Repeat for GALWEM ratios
     inquire(file=agrmet_struc(n)%galwem_nrt_bias_ratio_file, &
          exist=file_exists)
     if (.not. file_exists) then
        write(LIS_logunit,*) '[ERR] Cannot find ', &
             trim(agrmet_struc(n)%galwem_nrt_bias_ratio_file)
        call LIS_endrun()
     end if
     call LIS_verify( &
          nf90_open(path=agrmet_struc(n)%galwem_nrt_bias_ratio_file, &
          mode=NF90_NOWRITE, ncid=ncid), &
          'Error in nf90_open in USAF_pcpBackBiasRatio_nrt')
     call LIS_verify(nf90_inq_varid(ncid, "biasRatio", ppt_ratio_Id), &
          'biasRatio field not found in bias file')
     call LIS_verify(nf90_get_var(ncid, ppt_ratio_Id, &
          agrmet_struc(n)%galwem_nrt_bias_ratio, &
          start=(/1,1,imonth/), &
          count=(/LIS_rc%gnc(n),LIS_rc%gnr(n),1/)),&
          'n90_get_var failed for biasRatio in LIS param file')
     write(LIS_logunit,*) &
          '[INFO] Read GALWEM bias ratios for month ', imonth, &
          ' from ', trim(agrmet_struc(n)%galwem_nrt_bias_ratio_file)
     call LIS_verify(nf90_close(ncid),&
          'nf90_close failed in USAF_pcpBackBiasRatio_nrt')

#endif

   end subroutine USAF_pcpBackBiasRatio_nrt

   subroutine handle_mpi_error(errorcode, msg)
     use LIS_logMod, only:  LIS_logunit, LIS_endrun
     use LIS_mpiMod
     implicit none
     integer, intent(in) :: errorcode
     character(*), intent(in) :: msg
     character*(MPI_MAX_ERROR_STRING) :: buf
     integer :: resultlen, ierr
     if (errorcode .ne. MPI_SUCCESS) then
        write(LIS_logunit,*)'[ERR] ', trim(msg)
        call MPI_error_string(errorcode, buf, resultlen, ierr)
        write(LIS_logunit,*)'[ERR] MPI error: ', trim(buf)
        flush(LIS_logunit)
        call LIS_endrun()
     end if
   end subroutine handle_mpi_error
 end module USAF_bratsethMod
