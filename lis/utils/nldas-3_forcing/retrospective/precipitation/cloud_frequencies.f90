!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.8
!
! Copyright (c) 2026 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
!
! !ROUTINE: cloud_frequencies
!
! !REVISION HISTORY:
!  12 Jun 2026  Fadji Maina; Initial specification
!
! !COMPILATION:
!  This program can be compiled on NASA Discover using the Intel Fortran
!  compiler and the LISF NetCDF/HDF5 libraries as follows:
!
!  ifort -g -check all -traceback -names lowercase -convert big_endian 
!    -assume byterecl 
!    -I/discover/nobackup/projects/lis/libs/sles-12.3/netcdf/4.8.1_intel-2021.4.0/include 
!    cloud_frequencies.f90 -o cloud_frequencies 
!    -L/discover/nobackup/projects/lis/libs/sles-12.3/netcdf/4.8.1_intel-2021.4.0/lib 
!    -L/discover/nobackup/projects/lis/libs/sles-12.3/hdf5/1.12.1_intel-2021.4.0/lib 
!    -lnetcdff -lnetcdf -lhdf5_hl -lhdf5 -lz -lcurl 
!    -Wl,--no-relax -shared-intel
!
! !DESCRIPTION:
!  This program computes monthly cloud frequency at 1-km resolution from
!  daily MODIS MOD09GA surface reflectance QA files. The MODIS state_1km_1
!  quality-assurance variable is read for each daily file, and the internal
!  cloud algorithm flag is extracted from bit 10. Land/water QA bits are used
!  to mask ocean pixels and coastline pixels located within a specified search
!  radius of ocean water. Daily cloud flags are accumulated over each month and
!  converted to monthly cloud frequency using only valid land pixels. Monthly
!  values are written to NetCDF files together with the number of valid days,
!  coastline/ocean mask information, and latitude/longitude coordinates.
!
!  The output cloud frequency ranges from 0 to 1, where 0 indicates that no
!  valid daily observations were cloudy during the month and 1 indicates that
!  all valid daily observations were cloudy. Grid cells with insufficient valid
!  observations are assigned the specified fill value.
!
!EOP

program cloud_frequencies
    use netcdf
    implicit none

    ! ============================================================
    ! USER SETTINGS
    ! Change only this block
    ! ============================================================
    integer, parameter :: start_year  = 2026
    integer, parameter :: start_month = 3
    integer, parameter :: start_day   = 1

    integer, parameter :: end_year    = 2026
    integer, parameter :: end_month   = 12
    integer, parameter :: end_day     = 31

    character(len=*), parameter :: input_root  = "./MOD09GA.061_NLDAS3"
    character(len=*), parameter :: output_root = "./Cloudcover_1km"

    integer, parameter :: radius_pixels = 10

    ! ============================================================
    ! GRID SETTINGS
    ! ============================================================
    integer, parameter :: NX = 14040
    integer, parameter :: NY = 7800

    real, parameter :: lon_min = -169.0
    real, parameter :: lat_min = 7.0
    real, parameter :: dx = 0.00833
    real, parameter :: dy = 0.00833

    real, parameter :: fill_value = -9999.0

    ! MODIS QA bit settings
    integer, parameter :: cloud_bit_shift = 10
    integer, parameter :: landwater_bit_shift = 3
    integer, parameter :: landwater_mask = 7

    ! Minimum number of valid days required for monthly mean
    integer, parameter :: min_valid_days = 20

    ! ============================================================
    ! VARIABLES
    ! ============================================================
    integer :: i, j, ii, jj
    integer :: iyear, imonth, iday
    integer :: first_day, last_day
    integer :: doy
    integer :: dbin
    integer :: cloud_flag, lw_flag
    logical :: has_ocean

    character(len=256) :: input_file
    character(len=256) :: output_file
    character(len=16)  :: year_string
    character(len=16)  :: month_string
    character(len=16)  :: doy_string

    ! NetCDF variables
    integer :: ncid
    integer :: varid_state
    integer :: lat_dimid, lon_dimid
    integer :: lat_id, lon_id
    integer :: cfm_id, days_id, flag_id
    integer :: latitude_id, longitude_id
    integer :: dimids_2d(2)

    ! Arrays
    real, allocatable :: latitude(:,:), longitude(:,:)
    real, allocatable :: lat(:), lon(:)
    real, allocatable :: cfm_month(:,:), cfm_daily(:,:)
    real, allocatable :: count_day(:,:)
    real, allocatable :: coast_flag(:,:)

    integer, allocatable :: cfm_bin(:,:)
    integer, allocatable :: cell_day(:,:)

    ! ============================================================
    ! EXECUTABLE SECTION STARTS HERE
    ! Nothing should be declared below this point
    ! ============================================================

    allocate(latitude(NX,NY), longitude(NX,NY))
    allocate(lat(NY), lon(NX))

    allocate(cfm_bin(NX,NY))
    allocate(cfm_daily(NX,NY))
    allocate(cfm_month(NX,NY))
    allocate(cell_day(NX,NY))
    allocate(count_day(NX,NY))
    allocate(coast_flag(NX,NY))

    ! ============================================================
    ! BUILD LAT/LON ARRAYS
    ! ============================================================
    do i = 1, NX
        lon(i) = lon_min + dx * real(i - 1)
    end do

    do j = 1, NY
        lat(j) = lat_min + dy * real(j - 1)
    end do

    do i = 1, NX
        do j = 1, NY
            longitude(i,j) = lon(i)
            latitude(i,j)  = lat(j)
        end do
    end do

    ! ============================================================
    ! MAIN LOOP OVER YEARS AND MONTHS
    ! ============================================================
    do iyear = start_year, end_year

        do imonth = 1, 12

            if (.not. month_in_range(iyear, imonth, start_year, start_month, end_year, end_month)) cycle

            first_day = 1
            last_day  = days_in_month(iyear, imonth)

            if (iyear == start_year .and. imonth == start_month) first_day = start_day
            if (iyear == end_year   .and. imonth == end_month)   last_day  = end_day

            write(*,*) "Processing year/month: ", iyear, imonth

            cfm_month = 0.0
            cell_day  = 0

            do iday = first_day, last_day

                doy = day_of_year(iyear, imonth, iday)

                write(year_string, '(I4.4)') iyear
                write(month_string,'(I2.2)') imonth
                write(doy_string,  '(I3.3)') doy

                input_file = trim(input_root)//"/"//trim(year_string)// &
                             "/MOD09GA.061_"//trim(year_string)//trim(doy_string)//".nc4"

                write(*,*) "Reading: ", trim(input_file)

                ! ====================================================
                ! READ MODIS QA VARIABLE
                ! ====================================================
                call check(nf90_open(trim(input_file), NF90_NOWRITE, ncid))
                call check(nf90_inq_varid(ncid, "state_1km_1", varid_state))
                call check(nf90_get_var(ncid, varid_state, cfm_bin))
                call check(nf90_close(ncid))

                ! ====================================================
                ! STEP 1: EXTRACT CLOUD FLAG AND LAND/WATER FLAG
                ! ====================================================
                do i = 1, NX
                    do j = 1, NY

                        dbin = cfm_bin(i,j)

                        if (dbin == 65535) then
                            cfm_daily(i,j) = fill_value
                            coast_flag(i,j) = fill_value
                        else
                            ! Internal cloud algorithm flag: bit 10
                            cloud_flag = iand(ishft(dbin, -cloud_bit_shift), 1)
                            cfm_daily(i,j) = real(cloud_flag)

                            ! Land/water flag: bits 3-5
                            lw_flag = iand(ishft(dbin, -landwater_bit_shift), landwater_mask)

                            if (lw_flag == 2) then
                                coast_flag(i,j) = 1.0
                            else
                                coast_flag(i,j) = 0.0
                            end if
                        end if

                    end do
                end do

                ! ====================================================
                ! STEP 2: MASK OCEAN AND COASTLINE NEAR OCEAN
                ! ====================================================
                do i = 1, NX
                    do j = 1, NY

                        if (coast_flag(i,j) == 1.0) then

                            has_ocean = .false.

                            do ii = max(1, i - radius_pixels), min(NX, i + radius_pixels)
                                do jj = max(1, j - radius_pixels), min(NY, j + radius_pixels)

                                    lw_flag = iand(ishft(cfm_bin(ii,jj), -landwater_bit_shift), landwater_mask)

                                    if (lw_flag == 0 .or. lw_flag == 6 .or. lw_flag == 7) then
                                        has_ocean = .true.
                                    end if

                                end do
                            end do

                            if (has_ocean) then
                                cfm_daily(i,j) = fill_value
                            end if

                        end if

                        lw_flag = iand(ishft(cfm_bin(i,j), -landwater_bit_shift), landwater_mask)

                        if (lw_flag == 0 .or. lw_flag == 6 .or. lw_flag == 7) then
                            cfm_daily(i,j) = fill_value
                        end if

                    end do
                end do

                ! ====================================================
                ! STEP 3: ACCUMULATE DAILY CLOUD FLAGS INTO MONTHLY SUM
                ! ====================================================
                do i = 1, NX
                    do j = 1, NY

                        if (cfm_daily(i,j) >= 0.0) then
                            cfm_month(i,j) = cfm_month(i,j) + cfm_daily(i,j)
                            cell_day(i,j)  = cell_day(i,j) + 1
                        end if

                    end do
                end do

            end do   ! day loop

            ! ========================================================
            ! STEP 4: COMPUTE MONTHLY CLOUD FREQUENCY
            ! ========================================================
            do i = 1, NX
                do j = 1, NY

                    count_day(i,j) = real(cell_day(i,j))

                    if (cell_day(i,j) > min_valid_days .and. cfm_month(i,j) >= 0.0) then
                        cfm_month(i,j) = cfm_month(i,j) / real(cell_day(i,j))
                    else
                        cfm_month(i,j) = fill_value
                    end if

                    if (cell_day(i,j) == 0) then
                        cfm_month(i,j) = fill_value
                    end if

                end do
            end do

            ! ========================================================
            ! STEP 5: WRITE MONTHLY NETCDF OUTPUT
            ! ========================================================
            write(year_string, '(I4.4)') iyear
            write(month_string,'(I2.2)') imonth

            output_file = trim(output_root)//"/CFM_"//trim(year_string)//trim(month_string)//".nc"

            write(*,*) "Writing: ", trim(output_file)

            call check(nf90_create(trim(output_file), NF90_CLOBBER, ncid))

            call check(nf90_def_dim(ncid, "latitude",  NY, lat_dimid))
            call check(nf90_def_dim(ncid, "longitude", NX, lon_dimid))

            call check(nf90_def_var(ncid, "latitude",  NF90_REAL, lat_dimid, lat_id))
            call check(nf90_def_var(ncid, "longitude", NF90_REAL, lon_dimid, lon_id))

            dimids_2d = (/lon_dimid, lat_dimid/)

            call check(nf90_def_var(ncid, "CFM",  NF90_REAL, dimids_2d, cfm_id))
            call check(nf90_put_att(ncid, cfm_id, "_FillValue", fill_value))

            call check(nf90_def_var(ncid, "Days", NF90_REAL, dimids_2d, days_id))
            call check(nf90_put_att(ncid, days_id, "_FillValue", fill_value))

            call check(nf90_def_var(ncid, "Flag", NF90_REAL, dimids_2d, flag_id))
            call check(nf90_put_att(ncid, flag_id, "_FillValue", fill_value))

            call check(nf90_def_var(ncid, "lat", NF90_REAL, dimids_2d, latitude_id))
            call check(nf90_def_var(ncid, "lon", NF90_REAL, dimids_2d, longitude_id))

            call check(nf90_enddef(ncid))

            call check(nf90_put_var(ncid, lat_id, lat))
            call check(nf90_put_var(ncid, lon_id, lon))
            call check(nf90_put_var(ncid, cfm_id, cfm_month))
            call check(nf90_put_var(ncid, days_id, count_day))
            call check(nf90_put_var(ncid, flag_id, coast_flag))
            call check(nf90_put_var(ncid, latitude_id, latitude))
            call check(nf90_put_var(ncid, longitude_id, longitude))

            call check(nf90_close(ncid))

        end do   ! month loop

    end do   ! year loop

contains

    subroutine check(status)
        integer, intent(in) :: status

        if (status /= nf90_noerr) then
            print *, trim(nf90_strerror(status))
            stop "Stopped"
        end if
    end subroutine check


    logical function is_leap_year(year)
        integer, intent(in) :: year

        is_leap_year = .false.

        if (mod(year, 400) == 0) then
            is_leap_year = .true.
        else if (mod(year, 100) == 0) then
            is_leap_year = .false.
        else if (mod(year, 4) == 0) then
            is_leap_year = .true.
        end if

    end function is_leap_year


    integer function days_in_month(year, month)
        integer, intent(in) :: year, month

        select case(month)
        case(1,3,5,7,8,10,12)
            days_in_month = 31
        case(4,6,9,11)
            days_in_month = 30
        case(2)
            if (is_leap_year(year)) then
                days_in_month = 29
            else
                days_in_month = 28
            end if
        case default
            print *, "Invalid month: ", month
            stop
        end select

    end function days_in_month


    integer function day_of_year(year, month, day)
        integer, intent(in) :: year, month, day
        integer :: m

        day_of_year = day

        do m = 1, month - 1
            day_of_year = day_of_year + days_in_month(year, m)
        end do

    end function day_of_year


    logical function month_in_range(year, month, sy, sm, ey, em)
        integer, intent(in) :: year, month
        integer, intent(in) :: sy, sm, ey, em

        month_in_range = .true.

        if (year < sy) month_in_range = .false.
        if (year > ey) month_in_range = .false.

        if (year == sy .and. month < sm) month_in_range = .false.
        if (year == ey .and. month > em) month_in_range = .false.

    end function month_in_range

end program cloud_frequencies
