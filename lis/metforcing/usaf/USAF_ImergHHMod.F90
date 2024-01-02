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

! Macros for tracing - Requires ESMF 7_1_0+
#ifdef ESMF_TRACE
#define TRACE_ENTER(region) call ESMF_TraceRegionEnter(region)
#define TRACE_EXIT(region) call ESMF_TraceRegionExit(region)
#else
#define TRACE_ENTER(region)
#define TRACE_EXIT(region)
#endif

! Module for processing NASA half-hourly IMERG 30-min precipitation data.
! Updated 26 Apr 2022 by Eric Kemp/SSAI, to reduce memory footprint.
! Updated 14 Jul 2022 by Eric Kemp/SSAI, to support IMERG V07
! Updated 24 Aug 2023 by Eric Kemp/SSAI, to add alert files if problem
!   occurs with fetching IMERG.
module USAF_ImergHHMod

   ! Imports
   use ESMF

   ! Defaults
   implicit none
   private

   ! Public type for storing 30-min data and converting to 3-hr accums.
   type ImergHHPrecip
      private
      integer :: nlons
      integer :: nlats
      integer :: ntimes
      real, allocatable :: precip_cal_3hr(:,:)     ! mm
      real :: swlat
      real :: swlon
      real :: dlat
      real :: dlon
   end type ImergHHPrecip
   public :: ImergHHPrecip

   ! Public methods
   public :: newImergHHPrecip
   public :: destroyImergHHPrecip
   public :: update30minImergHHPrecip
   public :: count3hrObsImergHHPrecip
   public :: copyToObsDataImergHHPrecip
   public :: create_Imerg_HH_filename
   public :: fetch3hrImergHH

contains

   ! Constructor
   function newImergHHPrecip() result(this)
      implicit none
      type(ImergHHPrecip) :: this
      integer :: nlats, nlons, ntimes
      TRACE_ENTER("newImergHHPrecip")
      nlats = 1800
      nlons = 3600
      ntimes = 6 ! 6 30-min periods
      this%nlats = nlats
      this%nlons = nlons
      this%ntimes = ntimes
      !NOTE:  IMERG HDF5 grids are dimensioned lat,lon
      allocate(this%precip_cal_3hr(nlats, nlons))
      this%precip_cal_3hr(:,:) = 0
      this%swlat =  -89.95
      this%swlon = -179.95
      this%dlat = 0.1
      this%dlon = 0.1
      TRACE_EXIT("newImergHHPrecip")
   end function newImergHHPrecip

   ! Destructor
   subroutine destroyImergHHPrecip(this)
      implicit none
      type(ImergHHPrecip), intent(inout) :: this
      TRACE_ENTER("destroyImergHHPrecip")
      this%nlats = 0
      this%nlons = 0
      this%ntimes = 0
      this%swlat = 0
      this%swlon = 0
      this%dlat = 0
      this%dlon = 0
      if (allocated(this%precip_cal_3hr)) &
           deallocate(this%precip_cal_3hr)
      TRACE_EXIT("destroyImergHHPrecip")
   end subroutine destroyImergHHPrecip

   ! Copy to obsData
   subroutine copyToObsDataImergHHPrecip(this, sigmaOSqr, &
        oErrScaleLength, net, platform, obsData_struc)

      ! Modules
      use USAF_bratsethMod, only: USAF_obsData, USAF_assignObsData

      ! Defaults
      implicit none

      ! Arguments
      type(ImergHHPrecip), intent(in) :: this
      real, intent(in) :: sigmaOSqr
      real, intent(in) :: oErrScaleLength
      character*32, intent(in) :: net
      character*32, intent(in) :: platform
      type(USAF_ObsData), intent(inout) :: obsData_struc

      ! Local variables
      real :: ob, lat, lon
      integer :: i, j

      TRACE_ENTER("copyToObsDataImergHHPrecip")

      do j = 1, this%nlons
         do i = 1, this%nlats
            ob = this%precip_cal_3hr(i,j)
            if (ob .lt. 0) cycle

            lat = (this%swlat) + (i-1)*(this%dlat)
            lon = (this%swlon) + (j-1)*(this%dlon)
            call USAF_assignObsData(obsData_struc, net, platform,&
                 ob, lat, lon, sigmaOSqr, oErrScaleLength)

         end do ! i
      end do ! j

      TRACE_EXIT("copyToObsDataImergHHPrecip")
   end subroutine copyToObsDataImergHHPrecip

   ! Read slice from IMERG-E 30-min file.  Most low-level work occurs
   ! in internal subroutines.
   ! Code is designed to allow LIS to gracefully handle problems with
   ! HDF5 file.
   subroutine update30minImergHHPrecip(this, itime, filename, &
        plp_thresh, version)

      ! Imports
#if (defined USE_HDF5)
      use HDF5
#endif
      use LIS_coreMod, only: LIS_masterproc
      use LIS_logMod, only:  LIS_logunit, LIS_abort, LIS_endrun, &
           LIS_alert
      use LIS_mpiMod

      ! Defaults
      implicit none

      ! Arguments
      type(ImergHHPrecip), intent(inout) :: this
      integer, intent(in) :: itime
      character(len=*), intent(in) :: filename
      integer*2, intent(in) :: plp_thresh
      character(*), intent(in) :: version

      ! Local variables
      logical :: fail
      integer :: hdferr
#if (defined USE_HDF5)
      integer(HID_T) :: file_id, dataset_id, datatype_id
      integer(HSIZE_T) :: dims(3)
#endif
      integer :: i, j
      real, allocatable :: tmp_precip_cal(:,:,:)
      integer*2, allocatable :: tmp_prob_liq_precip(:,:,:)
      integer*2, allocatable :: tmp_ir_kalman_weights(:,:,:)
      integer :: icount
      character(len=255) :: message(20)
      integer :: ierr
      logical :: saved_good
      logical :: version_good
      character(22) :: varname
      logical :: apply_irkalman_test
      integer, save :: alert_number = 1
      logical :: file_exists

      TRACE_ENTER("update30minImergHHPrecip")

      saved_good = .false. ! Updated below

!Only define actual subroutine if LIS was compiled with HDF5 support.
#if (defined USE_HDF5)

      ! Only IMERG V06 and V07 supported
      version_good = .false.
      if (index(trim(version), 'V06') .ne. 0) then
         version_good = .true.
      else if (index(trim(version), 'V07') .ne. 0) then
         version_good = .true.
      end if
      if (.not. version_good) then
         write(LIS_logunit,*)&
              '[ERR] update30minImergHHPrecip Invalid IMERG Version  ', &
              trim(version)
         write(LIS_logunit,*) &
              'Only version generations 6 and 7 supported! '
         flush(LIS_logunit)
         message(:) = ''
         message(1) = '[ERR] Program:  LIS'
         message(2) = '  Routine: update30minImergHHPrecip.'
         message(3) = '  Invalid IMERG Version ' // trim(version)
         if(LIS_masterproc) then
            call LIS_alert( 'LIS.update30minImergHHPrecip', &
                 alert_number, &
                 message )
            alert_number = alert_number + 1
            call LIS_abort( message)
         endif
#if (defined SPMD)
         call MPI_Barrier(LIS_MPI_COMM, ierr)
#endif
      end if

      ! Sanity checks
      if (itime .lt. 1 .or. itime .gt. this%ntimes) then
         write(LIS_logunit,*)&
              '[ERR] update30minImergHHPrecip Invalid time level ', itime
         write(LIS_logunit,*) 'Must be in range from 1 to ', this%ntimes
         write(LIS_logunit,*)'Must be in range from 1 to ', this%ntimes
         flush(LIS_logunit)
         message(:) = ''
         message(1) = '[ERR] Program:  LIS'
         message(2) = '  Routine: update30minImergHHPrecip.'
         message(3) = '  Invalid time level selected'
         if(LIS_masterproc) then
            call LIS_alert( 'LIS.update30minImergHHPrecip', &
                 alert_number, &
                 message )
            alert_number = alert_number + 1
            call LIS_abort( message)
         endif
#if (defined SPMD)
         call MPI_Barrier(LIS_MPI_COMM, ierr)
#endif
      end if

      ! Initialize IDs.  Useful later for error handling.
      file_id = -1
      dataset_id = -1
      datatype_id = -1

      ! Initialize HDF5 Fortran interface.
      call open_hdf5_f_interface(fail)
      if (fail) then
         write(LIS_logunit,*)&
              '[WARN] update30minImergHHPrecip Cannot open HDF5 interface'
         flush(LIS_logunit)
         message(:) = ''
         message(1) = '[ERR] Program:  LIS'
         message(2) = '  Routine: update30minImergHHPrecip.'
         message(3) = '  Cannot open HDF5 interface'
         if(LIS_masterproc) then
            call LIS_alert( 'LIS.update30minImergHHPrecip', &
                 alert_number, &
                 message )
         endif
         alert_number = alert_number + 1
         goto 100
      end if

      ! Open the file
      call open_imerg_file(trim(filename), file_id, fail)
      if (fail) then
         write(LIS_logunit,*)&
              '[WARN] update30minImergHHPrecip Cannot open IMERG file', &
              trim(filename)
         flush(LIS_logunit)
         message(:) = ''
         message(1) = '[ERR] Program:  LIS'
         message(2) = '  Routine: update30minImergHHPrecip.'
         message(3) = '  Cannot open IMERG file ' // trim(filename)
         if(LIS_masterproc) then
            call LIS_alert( 'LIS.update30minImergHHPrecip', &
                 alert_number, &
                 message )
         endif
         alert_number = alert_number + 1
         goto 100
      end if

      ! Open the precipitationCal dataset; sanity check the data type,
      ! dimensions, and units; then read it in.
      ! EMK 14 Jul 2022 -- Support IMERG V06 or V07.
      if (index(trim(version), "V06") .ne. 0) then
         varname = "/Grid/precipitationCal"
      else if (index(trim(version), "V07") .ne. 0) then
         varname = "/Grid/precipitation"
      end if
      call open_imerg_dataset(file_id, trim(varname), &
           dataset_id, fail)
      if (fail) then
         write(LIS_logunit,*)&
              '[WARN] update30minImergHHPrecip Cannot open dataset ', &
              'in ', trim(filename)
         flush(LIS_logunit)
         message(:) = ''
         message(1) = '[ERR] Program:  LIS'
         message(2) = '  Routine: update30minImergHHPrecip.'
         message(3) = '  Cannot open dataset in ' // trim(filename)
         if(LIS_masterproc) then
            call LIS_alert( 'LIS.update30minImergHHPrecip', &
                 alert_number, &
                 message )
         endif
         alert_number = alert_number + 1
         goto 100
      end if

      call get_imerg_datatype(dataset_id, datatype_id, fail)
      if (fail) then
         write(LIS_logunit,*)&
              '[WARN] update30minImergHHPrecip Cannot get datatype ', &
              'in ', trim(filename)
         flush(LIS_logunit)
         message(:) = ''
         message(1) = '[ERR] Program:  LIS'
         message(2) = '  Routine: update30minImergHHPrecip.'
         message(3) = '  Cannot get datatype in ' // trim(filename)
         if(LIS_masterproc) then
            call LIS_alert( 'LIS.update30minImergHHPrecip', &
                 alert_number, &
                 message )
         endif
         alert_number = alert_number + 1
         goto 100
      end if

      call check_imerg_type(datatype_id, H5T_IEEE_F32LE, fail)
      if (fail) then
         write(LIS_logunit,*)&
              '[WARN] update30minImergHHPrecip Found bad datatype ', &
              'in ', trim(filename)
         flush(LIS_logunit)
         message(:) = ''
         message(1) = '[ERR] Program:  LIS'
         message(2) = '  Routine: update30minImergHHPrecip.'
         message(3) = '  Found bad datatype in ' // trim(filename)
         if(LIS_masterproc) then
            call LIS_alert( 'LIS.update30minImergHHPrecip', &
                 alert_number, &
                 message )
         endif
         alert_number = alert_number + 1
         goto 100
      end if

      dims(1) = this%nlats
      dims(2) = this%nlons
      dims(3) = 1
      call check_imerg_dims(dataset_id, 3, dims, fail)
      if (fail) then
         write(LIS_logunit,*)&
              '[WARN] update30minImergHHPrecip Found bad dimension ', &
              'in ', trim(filename)
         flush(LIS_logunit)
         message(:) = ''
         message(1) = '[ERR] Program:  LIS'
         message(2) = '  Routine: update30minImergHHPrecip.'
         message(3) = '  Found bad dimension in ' // trim(filename)
         if(LIS_masterproc) then
            call LIS_alert( 'LIS.update30minImergHHPrecip', &
                 alert_number, &
                 message )
         endif
         alert_number = alert_number + 1
         goto 100
      end if

      call check_imerg_units(dataset_id, "mm/hr", fail)
      if (fail) then
         write(LIS_logunit,*)&
              '[WARN] update30minImergHHPrecip Found bad units ', &
              'in ', trim(filename)
         flush(LIS_logunit)
         message(:) = ''
         message(1) = '[ERR] Program:  LIS'
         message(2) = '  Routine: update30minImergHHPrecip.'
         message(3) = '  Found bad units in ' // trim(filename)
         if(LIS_masterproc) then
            call LIS_alert( 'LIS.update30minImergHHPrecip', &
                 alert_number, &
                 message )
         endif
         alert_number = alert_number + 1
         goto 100
      end if

      allocate(tmp_precip_cal(dims(1), dims(2), dims(3)))
      tmp_precip_cal = 0
      call h5dread_f(dataset_id, H5T_IEEE_F32LE, tmp_precip_cal, dims, &
           hdferr)
      if (hdferr .ne. 0) then
         write(LIS_logunit,*)&
              '[WARN] update30minImergHHPrecip Bad data read from ', &
              'in ', trim(filename)
         flush(LIS_logunit)
         message(:) = ''
         message(1) = '[ERR] Program:  LIS'
         message(2) = '  Routine: update30minImergHHPrecip.'
         message(3) = '  Bad data read from ' // trim(filename)
         if(LIS_masterproc) then
            call LIS_alert( 'LIS.update30minImergHHPrecip', &
                 alert_number, &
                 message )
         endif
         alert_number = alert_number + 1
         goto 100
      end if

      ! Close the precipitationCal types.
      if (datatype_id .gt. -1) call close_imerg_datatype(datatype_id, fail)
      if (dataset_id .gt. -1) call close_imerg_dataset(dataset_id, fail)

      ! Open the probabilityLiquidPrecipitation dataset; sanity check
      ! the data type and dimensions; then read it in.
      call open_imerg_dataset(file_id, &
           "/Grid/probabilityLiquidPrecipitation", &
           dataset_id, fail)
      if (fail) then
         write(LIS_logunit,*)&
              '[WARN] update30minImergHHPrecip Cannot open dataset ', &
              'in ', trim(filename)
         flush(LIS_logunit)
         message(:) = ''
         message(1) = '[ERR] Program:  LIS'
         message(2) = '  Routine: update30minImergHHPrecip.'
         message(3) = '  Cannot open dataset in ' // trim(filename)
         if(LIS_masterproc) then
            call LIS_alert( 'LIS.update30minImergHHPrecip', &
                 alert_number, &
                 message )
         endif
         alert_number = alert_number + 1
         goto 100
      end if

      call get_imerg_datatype(dataset_id, datatype_id, fail)
      if (fail) then
         write(LIS_logunit,*)&
              '[WARN] update30minImergHHPrecip Cannot get datetype ', &
              'in ', trim(filename)
         flush(LIS_logunit)
         message(:) = ''
         message(1) = '[ERR] Program:  LIS'
         message(2) = '  Routine: update30minImergHHPrecip.'
         message(3) = '  Cannot get datatype in ' // trim(filename)
         if(LIS_masterproc) then
            call LIS_alert( 'LIS.update30minImergHHPrecip', &
                 alert_number, &
                 message )
         endif
         alert_number = alert_number + 1
         goto 100
      end if

      call check_imerg_type(datatype_id, H5T_STD_I16LE, fail)
      if (fail) then
         write(LIS_logunit,*)&
              '[WARN] update30minImergHHPrecip Found bad datetype ', &
              'in ', trim(filename)
         flush(LIS_logunit)
         message(:) = ''
         message(1) = '[ERR] Program:  LIS'
         message(2) = '  Routine: update30minImergHHPrecip.'
         message(3) = '  Found bad datatype in ' // trim(filename)
         if(LIS_masterproc) then
            call LIS_alert( 'LIS.update30minImergHHPrecip', &
                 alert_number, &
                 message )
         endif
         alert_number = alert_number + 1
         goto 100
      end if

      dims(1) = this%nlats
      dims(2) = this%nlons
      dims(3) = 1
      call check_imerg_dims(dataset_id, 3, dims, fail)
      if (fail) then
         write(LIS_logunit,*)&
              '[WARN] update30minImergHHPrecip Found bad dimension ', &
              'in ', trim(filename)
         flush(LIS_logunit)
         message(:) = ''
         message(1) = '[ERR] Program:  LIS'
         message(2) = '  Routine: update30minImergHHPrecip.'
         message(3) = '  Found bad dimension in ' // trim(filename)
         if(LIS_masterproc) then
            call LIS_alert( 'LIS.update30minImergHHPrecip', &
                 alert_number, &
                 message )
         endif
         alert_number = alert_number + 1
         goto 100
      end if

      allocate(tmp_prob_liq_precip(dims(1), dims(2), dims(3)))
      tmp_prob_liq_precip = 0
      call h5dread_f(dataset_id, H5T_STD_I16LE, tmp_prob_liq_precip, &
           dims, hdferr)
      if (hdferr .ne. 0) then
         write(LIS_logunit,*)&
              '[WARN] update30minImergHHPrecip Bad data read from ', &
              trim(filename)
         flush(LIS_logunit)
         message(:) = ''
         message(1) = '[ERR] Program:  LIS'
         message(2) = '  Routine: update30minImergHHPrecip.'
         message(3) = '  Bad data read from ' // trim(filename)
         if(LIS_masterproc) then
            call LIS_alert( 'LIS.update30minImergHHPrecip', &
                 alert_number, &
                 message )
         endif
         alert_number = alert_number + 1
         goto 100
      end if

      ! Close the probabilityLiquidPrecipitation types.
      if (datatype_id .gt. -1) call close_imerg_datatype(datatype_id, &
           fail)
      if (dataset_id .gt. -1) call close_imerg_dataset(dataset_id, fail)

      ! Open the /Grid/IRkalmanFilterWeight dataset; sanity check the data
      ! type and dimensions; then read it in.
      ! EMK 14 Jul 2022 -- IRkalmanFilterWeight is removed in V07.
      apply_irkalman_test = .true.
      if (index(version, "V06") .ne. 0) then
         call open_imerg_dataset(file_id, "/Grid/IRkalmanFilterWeight", &
              dataset_id, fail)
         if (fail) then
            write(LIS_logunit,*)&
                 '[WARN] update30minImergHHPrecip Cannot open dataset ', &
                 'from ', trim(filename)
            flush(LIS_logunit)
            message(:) = ''
            message(1) = '[ERR] Program:  LIS'
            message(2) = '  Routine: update30minImergHHPrecip.'
            message(3) = '  Cannot open dataset from ' // trim(filename)
            if(LIS_masterproc) then
               call LIS_alert( 'LIS.update30minImergHHPrecip', &
                    alert_number, &
                    message )
            endif
            alert_number = alert_number + 1
            goto 100
         end if

         call get_imerg_datatype(dataset_id, datatype_id, fail)
         if (fail) then
            write(LIS_logunit,*)&
                 '[WARN] update30minImergHHPrecip Cannot get datatype ', &
                 'from ', trim(filename)
            flush(LIS_logunit)
            message(:) = ''
            message(1) = '[ERR] Program:  LIS'
            message(2) = '  Routine: update30minImergHHPrecip.'
            message(3) = '  Cannot get datatype from ' // trim(filename)
            if(LIS_masterproc) then
               call LIS_alert( 'LIS.update30minImergHHPrecip', &
                    alert_number, &
                    message )
            endif
            alert_number = alert_number + 1
            goto 100
         end if

         call check_imerg_type(datatype_id, H5T_STD_I16LE, fail)
         if (fail) then
            write(LIS_logunit,*)&
                 '[WARN] update30minImergHHPrecip Found bad datatype ', &
                 'in ', trim(filename)
            flush(LIS_logunit)
            message(:) = ''
            message(1) = '[ERR] Program:  LIS'
            message(2) = '  Routine: update30minImergHHPrecip.'
            message(3) = '  Found bad datatype in ' // trim(filename)
            if(LIS_masterproc) then
               call LIS_alert( 'LIS.update30minImergHHPrecip', &
                    alert_number, &
                    message )
            endif
            alert_number = alert_number + 1
            goto 100
         end if

         dims(1) = this%nlats
         dims(2) = this%nlons
         dims(3) = 1
         call check_imerg_dims(dataset_id, 3, dims, fail)
         if (fail) then
            write(LIS_logunit,*)&
                 '[WARN] update30minImergHHPrecip Found bad dimension ', &
                 'in ', trim(filename)
            flush(LIS_logunit)
            message(:) = ''
            message(1) = '[ERR] Program:  LIS'
            message(2) = '  Routine: update30minImergHHPrecip.'
            message(3) = '  Found bad dimension in ' // trim(filename)
            if(LIS_masterproc) then
               call LIS_alert( 'LIS.update30minImergHHPrecip', &
                    alert_number, &
                    message )
            endif
            alert_number = alert_number + 1
            goto 100
         end if

         allocate(tmp_ir_kalman_weights(dims(1), dims(2), dims(3)))
         tmp_ir_kalman_weights = 0
         call h5dread_f(dataset_id, H5T_STD_I16LE, &
              tmp_ir_kalman_weights, dims, hdferr)
         if (hdferr .ne. 0) then
            write(LIS_logunit,*) &
                 '[WARN] update30minImergHHPrecip Cannot read data ', &
                 'from ', trim(filename)
            flush(LIS_logunit)
            message(:) = ''
            message(1) = '[ERR] Program:  LIS'
            message(2) = '  Routine: update30minImergHHPrecip.'
            message(3) = '  Cannot read data from ' // trim(filename)
            if(LIS_masterproc) then
               call LIS_alert( 'LIS.update30minImergHHPrecip', &
                    alert_number, &
                    message )
            endif
            alert_number = alert_number + 1
            goto 100
         end if

         ! Close the IRkalmanFilterWeight types.
         if (datatype_id .gt. -1) call close_imerg_datatype(datatype_id, &
              fail)
         if (dataset_id .gt. -1) call close_imerg_dataset(dataset_id, &
              fail)
      else
         apply_irkalman_test = .false. ! EMK for IMERG V07
      end if

      ! Save the "good" precipitationCal data.
      ! Precipitation units are converted from rate (mm/hr) to
      ! accumulation (mm).
      saved_good = .true.
      icount = 0
      do j = 1,this%nlons
         do i = 1,this%nlats

            ! Reject if we are missing data at an earlier time.
            if (this%precip_cal_3hr(i,j) < 0) cycle

            ! Gross error checks
            if (tmp_precip_cal(i,j,1) < 0 .or. &
                 tmp_prob_liq_precip(i,j,1) < plp_thresh) then
               this%precip_cal_3hr(i,j) = -9999
               cycle
            end if

            ! EMK: IR Kalman Filter check for IMERG V06
            if (apply_irkalman_test) then
               if (tmp_ir_kalman_weights(i,j,1) < 0) then
                  this%precip_cal_3hr(i,j) = -9999
                  cycle
               end if
            end if

            ! Estimate is good.
            this%precip_cal_3hr(i,j) = this%precip_cal_3hr(i,j) + &
                 (tmp_precip_cal(i,j,1) * 0.5)

            icount = icount + 1
         end do ! i
      end do ! j

      write(LIS_logunit,*) &
           '[INFO] update30minImergHHPrecip found ', icount, &
           ' good 30-min calibrated estimates'

      ! Clean up before returning.
100   continue
      if (.not. saved_good) this%precip_cal_3hr = -9999
      if (allocated(tmp_precip_cal)) deallocate(tmp_precip_cal)
      if (allocated(tmp_prob_liq_precip)) deallocate(tmp_prob_liq_precip)
      if (allocated(tmp_ir_kalman_weights)) &
           deallocate(tmp_ir_kalman_weights)
      if (datatype_id .gt. -1) call close_imerg_datatype(datatype_id, &
           fail)
      if (dataset_id .gt. -1) call close_imerg_dataset(dataset_id, fail)
      if (file_id .gt. -1) call close_imerg_file(file_id, fail)
      call close_hdf5_f_interface(fail)

!If LIS was compiled without HDF5 support, have the subroutine print/log
!an error message and abort.
#else
      write(LIS_logunit,*) &
           '[ERR] Cannot read IMERG data unless LIS is compiled with HDF5!'
      write(LIS_logunit,*) &
           'Recompile with HDF5 and try again!'
      write(LIS_logunit,*)'ABORTING'
      flush(LIS_logunit)
      message(:) = ''
      message(1) = '[ERR] Program: LIS'
      message(2) = '  Routine update30minImergHHPrecip.'
      message(3) = '  LIS was not compiled with HDF5 support'
      if (LIS_masterproc) then
         call LIS_alert('LIS.update30minImergHHPrecip', alert_number, &
              message)
         alert_number = alert_number + 1
         call LIS_abort(message)
      end if
#if (defined SPMD)
      call MPI_Barrier(LIS_MPI_COMM, ierr)
#endif

#endif

      TRACE_EXIT("update30minImergHHPrecip")

   contains

#if (defined USE_HDF5)
      ! Internal subroutine.  Open the HDF5 Fortran interface
      subroutine open_hdf5_f_interface(fail)
         use HDF5
         use LIS_logMod, only: LIS_logunit
         implicit none
         logical, intent(out) :: fail
         integer :: hdferr
         fail = .false.
         call h5open_f(hdferr)
         if (hdferr .ne. 0) then
            write(LIS_logunit,*) &
                 '[ERR] update30minImergHHPrecip cannot initialize HDF5 ', &
                 'Fortran interface!'
            fail = .true.
         end if
      end subroutine open_hdf5_f_interface

      ! Internal subroutine.  Open the IMERG HDF5 file
      subroutine open_imerg_file(filename, file_id, fail)
         use HDF5
         use LIS_logMod, only: LIS_logunit
         implicit none
         character(len=*), intent(in) :: filename
         integer(HID_T), intent(out) :: file_id
         logical, intent(out) :: fail
         integer :: hdferr
         fail = .false.
         call h5fopen_f(trim(filename), H5F_ACC_RDONLY_F, file_id, hdferr)
         if (hdferr .ne. 0) then
            write(LIS_logunit,*) &
                 '[ERR] update30minImergHHPrecip cannot open file ', &
                 trim(filename)
            fail = .true.
         else
            write(LIS_logunit,*) &
                 '[INFO] Opened Imerg file ', trim(filename)
         end if
      end subroutine open_imerg_file

      ! Internal subroutine.  Open HDF5 dataset.
      subroutine open_imerg_dataset(file_id, dataset_name, dataset_id, &
           fail)
         use HDF5
         use LIS_logMod, only: LIS_logunit
         implicit none
         integer(HID_T), intent(in) :: file_id
         character(len=*), intent(in) :: dataset_name
         integer(HID_T), intent(out) :: dataset_id
         logical, intent(out) :: fail
         fail = .false.
         call h5dopen_f(file_id, trim(dataset_name), dataset_id, hdferr)
         if (hdferr .ne. 0) then
            write(LIS_logunit,*)&
                 '[ERR] update30minImergHHPrecip cannot open dataset ', &
                 trim(dataset_name)
            fail = .true.
         end if
      end subroutine open_imerg_dataset

      ! Internal subroutine.  Sanity check IMERG precipitation units.
      subroutine check_imerg_units(dataset_id, units, fail)

         ! Imports
         use HDF5
         use ISO_C_BINDING
         use LIS_logMod, only: LIS_logunit

         ! Defaults
         implicit none

         ! Arguments
         integer(HID_T), intent(in) :: dataset_id
         character(len=*), intent(in) :: units
         logical,intent(out) :: fail

         ! Local variables
         integer(HID_T) :: attr_id, type_id, space_id, memtype_id
         integer :: hdferr
         integer(size_t) :: size
         integer(SIZE_T), parameter :: sdim = 5
         integer(HSIZE_T), dimension(1:1) :: dims = (/1/)
         integer(HSIZE_T), dimension(1:1) :: maxdims
         character(len=sdim), dimension(:), allocatable, target :: rdata
         type(C_PTR) :: f_ptr
         integer :: i

         fail = .false.

         ! Open the attribute
         call h5aopen_f(dataset_id, 'Units', attr_id, hdferr)
         if (hdferr .ne. 0) then
            write(LIS_logunit,*) &
                 '[ERR] update30minImergHHPrecip cannot open attribute'
            fail = .true.
            return
         end if

         ! Get the attribute datatype
         call h5aget_type_f(attr_id, type_id, hdferr)
         if (hdferr .ne. 0) then
            write(LIS_logunit,*) &
                 '[ERR] update30minImergHHPrecip cannot get attribute datatype'
            call h5aclose_f(attr_id, hdferr)
            fail = .true.
            return
         end if

         ! Get the size of the attribute datatype, and sanity check.
         call h5tget_size_f(type_id, size, hdferr)
         if (hdferr .ne. 0) then
            write(LIS_logunit,*) &
                 '[ERR] update30minImergHHPrecip cannot get attribute ', &
                 'datatype size'
            call h5tclose_f(type_id, hdferr)
            call h5aclose_f(attr_id, hdferr)
            fail = .true.
            return
         end if
         if (size .gt. sdim+1) then
            write(LIS_logunit,*) &
                 '[ERR] update30minImergHHPrecip expected smaller attribute',&
                 'datatype size'
            write(LIS_logunit,*)'Expected ', sdim+1
            write(LIS_logunit,*)'Found ', size
            call h5tclose_f(type_id, hdferr)
            call h5aclose_f(attr_id, hdferr)
            fail = .true.
            return
         end if

         ! Get the attribute dataspace
         call h5aget_space_f(attr_id, space_id, hdferr)
         if (hdferr .ne. 0) then
            write(LIS_logunit,*) &
                 '[ERR] update30minImergHHPrecip cannot get attribute', &
                 'dataspace'
            call h5tclose_f(type_id, hdferr)
            call h5aclose_f(attr_id, hdferr)
            fail = .true.
            return
         end if

         ! Get the dimensions of the dataspace
         call h5sget_simple_extent_dims_f(space_id, dims, maxdims, hdferr)
         if (hdferr .ne. 0) then
            write(LIS_logunit,*) &
                 '[ERR] update30minImergHHPrecip cannot get attribute ', &
                 'dataspace dimensions'
            call h5sclose_f(space_id, hdferr)
            call h5tclose_f(type_id, hdferr)
            call h5aclose_f(attr_id, hdferr)
            fail = .true.
            return
         end if

         ! Create the memory datatype
         call h5tcopy_f(H5T_FORTRAN_S1, memtype_id, hdferr)
         if (hdferr .ne. 0) then
            write(LIS_logunit,*) &
                 '[ERR] update30minImergHHPrecip cannot copy attribute ', &
                 'memory datatype.'
            call h5sclose_f(space_id, hdferr)
            call h5tclose_f(type_id, hdferr)
            call h5aclose_f(attr_id, hdferr)
            fail = .true.
            return
         end if
         call h5tset_size_f(memtype_id, sdim, hdferr)
         if (hdferr .ne. 0) then
            write(LIS_logunit,*) &
                 '[ERR] update30minImergHHPrecip cannot set attribute ', &
                 'memory datatype size.'
            call h5tclose_f(memtype_id, hdferr)
            call h5sclose_f(space_id, hdferr)
            call h5tclose_f(type_id, hdferr)
            call h5aclose_f(attr_id, hdferr)
            fail = .true.
            return
         end if

         ! Read the attribute
         allocate(rdata(1:dims(1)))
         f_ptr = C_LOC(rdata(1)(1:1))
         call h5aread_f(attr_id, memtype_id, f_ptr, hdferr)
         if (hdferr .ne. 0) then
            write(LIS_logunit,*) &
                 '[ERR] update30minImergHHPrecip cannot read attribute.'
            deallocate(rdata)
            call h5tclose_f(memtype_id, hdferr)
            call h5sclose_f(space_id, hdferr)
            call h5tclose_f(type_id, hdferr)
            call h5aclose_f(attr_id, hdferr)
            fail = .true.
            return
         end if

         ! Check the units
         if (trim(rdata(1)) .ne. trim(units)) then
            write(LIS_logunit,*) &
                 '[ERR] update30minImergHHPrecip found wrong precipitation', &
                 'units'
            write(LIS_logunit,*) 'Expected mm/hr'
            write(LIS_logunit,*) 'Found ', trim(rdata(1))
            deallocate(rdata)
            call h5tclose_f(memtype_id, hdferr)
            call h5sclose_f(space_id, hdferr)
            call h5tclose_f(type_id, hdferr)
            call h5aclose_f(attr_id, hdferr)
            fail = .true.
            return
         end if

         ! Clean up
         deallocate(rdata)
         call h5tclose_f(memtype_id, hdferr)
         call h5sclose_f(space_id, hdferr)
         call h5tclose_f(type_id, hdferr)
         call h5aclose_f(attr_id, hdferr)

      end subroutine check_imerg_units

      ! Internal subroutine.  Get HDF5 datatype
      subroutine get_imerg_datatype(dataset_id, datatype_id, fail)
         use HDF5
         use LIS_logMod, only: LIS_logunit
         implicit none
         integer(HID_T), intent(in) :: dataset_id
         integer(HID_T), intent(out) :: datatype_id
         logical, intent(out) :: fail
         integer :: hdferr
         fail = .false.
         call h5dget_type_f(dataset_id, datatype_id, hdferr)
         if (hdferr .ne. 0) then
            write(LIS_logunit,*) &
                 '[ERR] update30minImergHHPrecip cannot determine datatype'
            fail = .true.
         end if
      end subroutine get_imerg_datatype

      ! Internal function.  Check datatype
      subroutine check_imerg_type(datatype_id, datatype, fail)
         use HDF5
         use LIS_logMod, only: LIS_logunit
         integer(HID_T), intent(in) :: datatype_id
         integer(HID_T), intent(in) :: datatype
         logical, intent(out) :: fail
         logical :: flag
         integer :: hdferr
         fail = .false.
         call h5tequal_f(datatype_id, datatype, flag, hdferr)
         if (hdferr .ne. 0) then
            write(LIS_logunit,*) &
                 '[ERR] update30minImergHHPrecip cannot confirm datatype!'
            fail = .true.
            return
         end if
         if (.not. flag) then
            write(LIS_logunit,*) &
                 '[ERR] update30minImergHHPrecip datatype is wrong type!'
            fail = .true.
            return
         end if
      end subroutine check_imerg_type

      ! Internal subroutine.  Check the rank/dimensions of dataset.
      subroutine check_imerg_dims(dataset_id, rank, dims, fail)

         ! Imports
         use HDF5
         use LIS_logMod, only: LIS_logunit

         ! Defaults
         implicit none

         ! Arguments
         integer(HID_T), intent(in) :: dataset_id
         integer, intent(in) :: rank
         integer(HSIZE_T), intent(in) :: dims(rank)
         logical, intent(out) :: fail

         ! Local variables
         integer(HID_T) :: dataspace_id
         integer :: hdferr, dataspace_rank
         integer(HSIZE_T), allocatable :: dataspace_dims(:)
         integer(HSIZE_T), allocatable :: dataspace_maxdims(:)
         logical :: flag
         integer :: i

         ! First, get the dataspace for the dataset
         call h5dget_space_f(dataset_id, dataspace_id, hdferr)
         if (hdferr .ne. 0) then
            write(LIS_logunit,*) &
                 '[ERR] update30minImergHHPrecip could not get dataspace'
            fail = .true.
            return
         end if

         ! Sanity check:  Make sure this dataspace is "simple"
         call h5sis_simple_f(dataspace_id, flag, hdferr)
         if (hdferr .ne. 0) then
            write(LIS_logunit,*) &
                 '[ERR] update30minImergHHPrecip cannot determine if ', &
                 'dataspace is simple'
            fail = .true.
            return
         end if
         if (.not. flag) then
            write(LIS_logunit,*) &
                 '[ERR] update30minImergHHPrecip dataspace is not simple'
            fail = .true.
            return
         end if

         ! Check the rank (number of dimensions)
         call h5sget_simple_extent_ndims_f(dataspace_id, dataspace_rank, &
              hdferr)
         if (hdferr .ne. 0) then
            write(LIS_logunit,*)&
                 '[ERR] update30minImergHHPrecip cannot get rank of dataspace '
            fail = .true.
            return
         end if
         if (dataspace_rank .ne. rank) then
            write(LIS_logunit,*) &
                 '[ERR] update30minImergHHPrecip expected rank ', rank
            write(LIS_logunit,*) 'But found rank ', dataspace_rank
            fail = .true.
            return
         end if

         ! Check the dimensions
         allocate(dataspace_dims(rank))
         allocate(dataspace_maxdims(rank))
         call h5sget_simple_extent_dims_f(dataspace_id, dataspace_dims, &
              dataspace_maxdims, hdferr)
         if (hdferr .ne. rank) then
            write(LIS_logunit,*) &
                 '[ERR] update30minImergHHPrecip cannot get dims for dataspace'
            deallocate(dataspace_dims)
            deallocate(dataspace_maxdims)
            fail = .true.
            return
         end if
         do i = 1, rank
            if (dataspace_dims(i) .ne. dims(i)) then
               fail = .true.
               exit
            end if
         end do
         if (fail) then
            write(LIS_logunit,*) &
              '[ERR] update30minImergHHPrecip found bad dimensions for ', &
              'dataspace'
            write(LIS_logunit,*) 'Expected ', dims(:)
            write(LIS_logunit,*) 'Found ', dataspace_dims(:)
            deallocate(dataspace_dims)
            deallocate(dataspace_maxdims)
            return
         end if

         ! Clean up
         deallocate(dataspace_dims)
         deallocate(dataspace_maxdims)
         call h5sclose_f(dataspace_id, hdferr)
         if (hdferr .ne. 0) then
            write(LIS_logunit,*) &
                 '[ERR] update30minImergHHPrecip cannot close dataspace'
            fail = .true.
            return
         end if

      end subroutine check_imerg_dims

      ! Internal subroutine.  Close datatype
      subroutine close_imerg_datatype(datatype_id, fail)
         use HDF5
         use LIS_logMod, only: LIS_logunit
         integer(HID_T), intent(inout) :: datatype_id
         logical, intent(out) :: fail
         integer :: hdferr
         fail = .false.
         call h5tclose_f(datatype_id, hdferr)
         if (hdferr .ne. 0) then
            write(LIS_logunit,*) &
                 '[ERR] update30minImergHHPrecip cannot close datatype '
            fail = .true.
         end if
         datatype_id = -1
      end subroutine close_imerg_datatype

      ! Internal function.  Close the dataset
      subroutine close_imerg_dataset(dataset_id, fail)
         use HDF5
         use LIS_logMod, only: LIS_logunit
         integer(HID_T), intent(inout) :: dataset_id
         logical, intent(out) :: fail
         integer :: hdferr
         fail = .false.
         call h5dclose_f(dataset_id, hdferr)
         if (hdferr .ne. 0) then
            write(LIS_logunit,*) &
                 '[ERR] update30minImergHHPrecip cannot close dataset '
            fail = .true.
         end if
         dataset_id = -1
      end subroutine close_imerg_dataset

      ! Internal subroutine.  Close the IMERG HDF5 file.
      subroutine close_imerg_file(file_id, fail)
         use HDF5
         use LIS_logMod, only: LIS_logunit
         implicit none
         integer(HID_T), intent(inout) :: file_id
         logical, intent(out) :: fail
         integer :: hdferr
         fail = .false.
         call h5fclose_f(file_id, hdferr)
         if (hdferr .ne. 0) then
            write(LIS_logunit,*) &
                 '[ERR] update30minImergHHPrecip cannot close file ', &
                 trim(filename)
            fail = .true.
         else
            write(LIS_logunit,*) &
                 '[INFO] Closed IMERG file ', trim(filename)
         end if
         file_id = -1
      end subroutine close_imerg_file

      ! Internal subroutine.  Close the HDF5 Fortran interface.
      subroutine close_hdf5_f_interface(fail)
         use HDF5
         use LIS_logMod, only: LIS_logunit
         implicit none
         logical, intent(out) :: fail
         integer :: hdferr
         fail = .false.
         call h5close_f(hdferr)
         if (hdferr .ne. 0) then
            write(LIS_logunit,*) &
                 '[ERR] update30minImergHHPrecip cannot close HDF5 Fortran ', &
                 'interface!'
            fail = .true.
         end if
      end subroutine close_hdf5_f_interface

#endif

   end subroutine update30minImergHHPrecip

   ! Determine number of Imerg points with valid data
   function count3hrObsImergHHPrecip(this) result(icount)

      ! Modules
      use LIS_logMod, only: LIS_logunit

      ! Defaults
      implicit none

      ! Arguments
      type(ImergHHPrecip), intent(in) :: this

      ! Result
      integer :: icount

      ! Local variables
      integer :: i,j

      TRACE_ENTER("count3hrObsImergHHPrecip")

      ! Pure calibrated PMW points.
      icount = 0
      do j = 1, this%nlons
         do i = 1, this%nlats
            if (this%precip_cal_3hr(i,j) .lt. 0) cycle
            icount = icount + 1
         end do ! i
      end do ! j
      write(LIS_logunit,*) &
           '[INFO] ImergHHPrecip found ', icount, &
           ' good 3-hr calibrated estimates.'

      TRACE_EXIT("count3hrObsImergHHPrecip")

   end function count3hrObsImergHHPrecip

   ! Construct IMERG 30-min HDF5 filename
   subroutine create_Imerg_HH_filename(dir, product, version,&
        yr, mo, da, hr, mn, filename)

      ! Imports
      use LIS_coreMod, only: LIS_masterproc
      use LIS_logMod, only:  LIS_logunit, LIS_abort, LIS_endrun, &
           LIS_alert
#if (defined SPMD)
      use LIS_mpiMod, only: LIS_MPI_COMM
#endif
      use LIS_timeMgrMod, only: LIS_calendar

      ! Defaults
      implicit none

      ! Arguments
      character(len=*), intent(in) :: dir
      character(len=*), intent(in) :: product
      character(len=*), intent(in) :: version
      integer, intent(in) :: yr
      integer, intent(in) :: mo
      integer, intent(in) :: da
      integer, intent(in) :: hr
      integer, intent(in) :: mn
      character(len=255), intent(out) :: filename

      ! Local variables
      integer :: tmp_yr, tmp_mo, tmp_da, tmp_hr, tmp_mn, tmp_ss
      integer :: tmp_minutes_in_day
      character(len=4) :: syr, sminutes_in_day
      character(len=2) :: smo, sda, shr, smn, sss
      type(ESMF_TIME) :: start_time, end_time, start_of_day
      type(ESMF_TIMEINTERVAL) :: half_hour
      type(ESMF_TIMEINTERVAL) :: time_diff
      character(len=255) :: message(20)
      integer :: ierr
      integer, save :: alert_number = 1
      logical :: file_exists
      integer :: rc

      TRACE_ENTER("create_Imerg_HH_filename")

      ! Construct start of filename.
      write(unit=syr, fmt='(i4.4)') yr
      write(unit=smo, fmt='(i2.2)') mo
      filename = trim(dir)//"/"//trim(syr)//trim(smo)//"/"
      filename = trim(filename)//trim(product)//".MS.MRG.3IMERG."

      ! Construct filename up through start date/time.
      write(unit=syr, fmt='(i4.4)') yr
      write(unit=smo, fmt='(i2.2)') mo
      write(unit=sda, fmt='(i2.2)') da
      write(unit=shr, fmt='(i2.2)') hr
      write(unit=smn, fmt='(i2.2)') mn
      write(unit=sss, fmt='(i2.2)') 0
      filename = trim(filename)//syr//smo//sda
      filename = trim(filename)//"-S"//shr//smn//sss

      ! Determine end date/time.
      call esmf_timeset(start_time, &
           yy=yr, &
           mm=mo, &
           dd=da, &
           h=hr, &
           m=mn, &
           s=0, &
           calendar = LIS_calendar, &
           rc = rc)
      call esmf_timeintervalset(half_hour, m=29, s=59, rc=rc)
      end_time = start_time + half_hour
      call esmf_timeget(end_time, &
           yy = tmp_yr, &
           mm = tmp_mo, &
           dd = tmp_da, &
           h  = tmp_hr, &
           m  = tmp_mn, &
           s  = tmp_ss, &
           rc = rc)

      ! Add end date/time to filename
      write(unit=syr, fmt='(i4.4)') tmp_yr
      write(unit=smo, fmt='(i2.2)') tmp_mo
      write(unit=sda, fmt='(i2.2)') tmp_da
      write(unit=shr, fmt='(i2.2)') tmp_hr
      write(unit=smn, fmt='(i2.2)') tmp_mn
      write(unit=sss, fmt='(i2.2)') tmp_ss
      filename = trim(filename)//"-E"//shr//smn//sss

      ! Determine number of minutes between start_of_day and start_time
      call esmf_timeset(start_of_day, &
           yy=yr, &
           mm=mo, &
           dd=da, &
           h=0, &
           m=0, &
           s=0, &
           calendar = LIS_calendar, &
           rc = rc)
      time_diff = start_time - start_of_day
      call esmf_timeintervalget(time_diff, m = tmp_minutes_in_day)

      ! Append minutes from start of day to filename
      write(unit=sminutes_in_day, fmt='(i4.4)') tmp_minutes_in_day
      filename = trim(filename)//"."//sminutes_in_day

      ! Finish filename construction
      ! EMK...Accomodate Final Run
      select case (trim(product))
      case ("3B-HHR")
         filename = trim(filename)//"."//trim(version)//".HDF5"
      case ("3B-HHR-E", "3B-HHR-L")
         filename = trim(filename)//"."//trim(version)//".RT-H5"
      case default
         write(LIS_logunit,*) &
              '[ERR] Invalid IMERG product!'
         write(LIS_logunit,*) &
              '[ERR] Valid options are 3B-HHR, 3B-HHR-E, 3B-HHR-L'
         write(LIS_logunit,*) &
              '[ERR] Found ',trim(product)
         flush(LIS_logunit)
         message(:) = ''
         message(1) = '[ERR] Program:  LIS'
         message(2) = '  Routine: create_imerg_HH_filename.'
         message(3) = '  Invalid IMERG product selected'

         if (LIS_masterproc) then
            call LIS_alert( 'LIS.create_imerg_HH_filename', &
                 alert_number, &
                 message )
            alert_number = alert_number + 1
            call LIS_abort( message)
         end if

#if (defined SPMD)
         call MPI_Barrier(LIS_MPI_COMM, ierr)
#endif

      end select

      TRACE_EXIT("create_Imerg_HH_filename")

   end subroutine create_Imerg_HH_filename

   ! Driver routine to fetch 3hr IMERG data for given start date.
   subroutine fetch3hrImergHH(j3hr, datadir, product, version, &
        plp_thresh, nest, sigmaOSqr, oErrScaleLength, net, platform, &
        precipObsData)

      ! Modules
      use LIS_coreMod, only: LIS_masterproc
      use LIS_logMod, only: LIS_logunit, LIS_alert
      use LIS_timeMgrMod, only: LIS_julhr_date, LIS_calendar
      use USAF_bratsethMod, only: USAF_ObsData, USAF_createObsData

      ! Defaults
      implicit none

      ! Arguments
      integer, intent(in) :: j3hr
      character(len=*),intent(in) :: datadir
      character(len=*),intent(in) :: product
      character(len=*),intent(in) :: version
      integer*2, intent(in) :: plp_thresh
      integer,intent(in) :: nest
      real, intent(in) :: sigmaOSqr
      real, intent(in) :: oErrScaleLength
      character(len=*),intent(in) :: net
      character(len=*),intent(in) :: platform
      type(USAF_ObsData), intent(out) :: precipObsData

      ! Local variables
      type(ESMF_TIME) :: start_time, cur_time
      type(ESMF_TIMEINTERVAL) :: half_hour
      type(ImergHHPrecip) :: imerg
      integer :: yr, mo, da, hr, mn
      integer :: itime
      character(len=255) :: filename
      character(255) :: message(20)
      integer, save :: alert_number = 1
      logical :: file_exists
      integer :: icount
      integer :: rc

      ! Save the start time
      call LIS_julhr_date(j3hr, yr, mo, da, hr)
      call esmf_timeset(start_time, &
           yy=yr, &
           mm=mo, &
           dd=da, &
           h=hr, &
           m=0, &
           s=0, &
           calendar = LIS_calendar, &
           rc = rc)
      call esmf_timeset(cur_time, &
           yy=yr, &
           mm=mo, &
           dd=da, &
           h=hr, &
           m=0, &
           s=0, &
           calendar = LIS_calendar, &
           rc = rc)

      ! Set the half hour time interval
      call esmf_timeintervalset(half_hour, m=30, rc=rc)

      ! Create the ImergHHPrecip object
      imerg = newImergHHPrecip()

      ! Loop through each 30-min period for 3-hr accumulations
      itime = 1
      do

         ! Get current filename
         call esmf_timeget(cur_time, &
              yy = yr, &
              mm = mo, &
              dd = da, &
              h  = hr, &
              m  = mn, &
              rc = rc)
         call create_Imerg_HH_filename(datadir, product, version, &
              yr,mo,da,hr,mn,filename)

         ! Report if the file doesn't exist
         inquire(file=trim(filename), exist=file_exists)
         if (.not. file_exists) then
            write(LIS_logunit,*)'[WARN] Cannot find IMERG file ', &
                 trim(filename)
            write(LIS_logunit,*) &
                 '[WARN] No IMERG data will be assimilated'
            flush(LIS_logunit)
            message(:) = ''
            message(1) = '[WARN] Program: LIS'
            message(2) = '  Routine fetch3hrImergHH.'
            message(3) = '  Cannot find IMERG file ' // trim(filename)
            if (LIS_masterproc) then
               call LIS_alert('LIS.fetch3hrImergHH', alert_number, &
                    message)
            end if
            alert_number = alert_number + 1
            imerg%precip_cal_3hr = -9999
         else
            ! Process the 30HH file
            call update30minImergHHPrecip(imerg, itime, filename, &
                 plp_thresh, version)
         end if

         ! Next cycle
         itime = itime + 1
         if (itime .gt. 6) exit
         cur_time = cur_time + half_hour

      end do

      ! Create obsData object.  For efficiency, allocate memory to match
      ! the total number of good 3-hr values
      icount = count3hrObsImergHHPrecip(imerg)
      call USAF_createObsData(precipObsData, nest, maxobs=icount)

      ! Copy the good values into the ObsData object
      call copyToObsDataImergHHPrecip(imerg, sigmaOSqr, &
           oErrScaleLength, &
           net, platform, precipObsData)

      ! Clean up
      call destroyImergHHPrecip(imerg)

   end subroutine fetch3hrImergHH
end module USAF_ImergHHMod
