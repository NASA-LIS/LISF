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
! !ROUTINE: semivar
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
!    semivar.f90 -o semivar 
!    -L/discover/nobackup/projects/lis/libs/sles-12.3/netcdf/4.8.1_intel-2021.4.0/lib 
!    -L/discover/nobackup/projects/lis/libs/sles-12.3/hdf5/1.12.1_intel-2021.4.0/lib 
!    -lnetcdff -lnetcdf -lhdf5_hl -lhdf5 -lz -lcurl 
!    -Wl,--no-relax -shared-intel
!
! !DESCRIPTION:
!  This program computes empirical semivariograms of daily precipitation
!  errors between gridded MERRA-2/LIS precipitation and station observations.
!  The output is a text file containing distance bins, semivariogram values,
!  and pair counts.
!
!EOP

program semivar
  use netcdf
  implicit none

  ! ---- Parameters ----
  integer, parameter :: NX=2926, NY=1626
  integer, parameter :: ndims=2, maxobs=50000, max_vario_bins=51
  integer, parameter :: grid_cell_km = 4
  integer, parameter :: max_distance_km = 500
  integer, parameter :: search_radius_cells = max_distance_km / grid_cell_km
  real(8), parameter :: spatial_tol = 0.02  ! ~ half a cell (2 km) in degrees

  ! ---- Variables ----
  integer :: i, j, i1, j1, iobs, ii1, ii2, jj1, jj2, index
  integer :: iyear, imonth, iday, nday(12), nvar_day, nlines
  integer :: bm1, bm2, nobs, nobsd, nobsd1, count_year,d1,jj
  integer :: ncid, totprecip_id,id,jd
  character(len=300) :: filename1a, filenameobs, fileplacea, fileplaceb,filenamevar
  real(8) :: latb, lonb, lata, lona, dist,vario_bin_dist
  integer :: actual_year

  ! ---- Arrays ----
  real(8), allocatable :: lat1(:,:), lon1(:,:), precepimerg(:,:), precep_back(:,:), precep_obs(:,:)
  real(8), allocatable :: diffij(:,:), vario_d(:,:), lat(:), lon(:), vario(:),diffij1(:,:)
  integer, allocatable :: icounts_vario_d(:,:), icounts_vario(:)
  real(8), allocatable :: lat_obs(:), lon_obs(:), precepobs(:)
  real(8), allocatable :: lat_obsd(:), lon_obsd(:), precepobsd(:)
  integer, allocatable :: stat(:), countd(:), dup(:), useo(:)
  character(len=300) :: yyyymm, yyyymmdd

  ! ---- Initialize ----
  vario_bin_dist = 10.0
  count_year = 0
  nvar_day = 0
  nday = (/31,28,31,30,31,30,31,31,30,31,30,31/)

  allocate(precepimerg(NX,NY), precep_back(NX,NY), precep_obs(NX,NY), diffij(NX,NY))
  allocate(lat1(NX,NY), lon1(NX,NY))
  allocate(lat_obs(maxobs), lon_obs(maxobs), precepobs(maxobs), diffij1(NX,NY))
  allocate(lat_obsd(maxobs), lon_obsd(maxobs), precepobsd(maxobs))
  allocate(stat(maxobs), countd(maxobs), dup(maxobs), useo(maxobs))
  allocate(vario_d(max_vario_bins,7000), icounts_vario_d(max_vario_bins,7000))
  allocate(vario(max_vario_bins), icounts_vario(max_vario_bins))

  do i = 1, NX
    do j = 1, NY
      lat1(i,j) = 7.0 + 0.04 * j
      lon1(i,j) = -169.0 + 0.04 * i
    end do
  end do

  do iyear = 2026, 2026
    actual_year = iyear

    if (mod(actual_year, 4) == 0 .and. (mod(actual_year, 100) /= 0 .or. mod(actual_year, 400) == 0)) then
      nday(2) = 29
    else
      nday(2) = 28
    end if

    write(fileplaceb, '(A,I4.4,A)') './obs_raw/download/', actual_year, '/'

        if (actual_year==2025) bm1=12
        if (actual_year==2026) bm1=2
    do imonth = bm1, 12
        write(*,*) iyear,imonth
    do iday = 1, nday(imonth)
        ! Format YYYYMM and YYYYMMDD
        write(yyyymm, '(I4.4,I2.2)') actual_year, imonth
        write(yyyymmdd, '(I4.4,I2.2,I2.2)') actual_year, imonth, iday

        ! Now construct full file path
        filename1a = './merra2_daily/SURFACEMODEL/' // trim(yyyymm) // '/LIS_HIST_' // trim(yyyymmdd) // '0000.d01.nc'

        ! Format the date
        write(yyyymmdd, '(I4.4,I2.2,I2.2)') actual_year, imonth, iday

        ! Concatenate to build filename
        filenameobs = trim(fileplaceb) // 'ncei_obs_' // trim(yyyymmdd) // '.txt'

        call check(nf90_open(filename1a, NF90_NOWRITE, ncid))
        call check(nf90_inq_varid(ncid, "TotalPrecip_tavg", totprecip_id))
        call check(nf90_get_var(ncid, totprecip_id, precepimerg))
        call check(nf90_close(ncid))

        call read_obs_data(trim(filenameobs), nlines, lat_obsd, lon_obsd, precepobsd, nobsd1)

        precep_back = precepimerg * 86400.0d0
        precep_obs = -9999.0d0
        useo = 0

              do iobs = 1, nobsd1
                id = int((lon_obsd(iobs) + 169.0) / 0.04)
                jd = int((lat_obsd(iobs) - 7.0) / 0.04)
                if (id >= 1 .and. id <= NX .and. jd >= 1 .and. jd <= NY) then
                if (precep_back(id,jd) >= 0.0d0) then
                    precep_obs(id,jd) = precepobsd(iobs)        !/10.0
                end if
                end if
              end do

              nvar_day = nvar_day + 1
              vario=0.0       ;       icounts_vario=0
! Semivariogram calculation
do i = 1, nx
  do j = 1, ny
    if (precep_back(i,j) >= 0.0d0 .and. precep_obs(i,j) >= 0.0d0) then
      diffij(i,j) = precep_obs(i,j) - precep_back(i,j)

      latb = lat1(i,j)
      lonb = lon1(i,j)
      ii1 = max(1, i - search_radius_cells)
      ii2 = min(nx, i + search_radius_cells)
      jj1 = max(1, j - search_radius_cells)
      jj2 = min(ny, j + search_radius_cells)

      do i1 = ii1, ii2
        do j1 = jj1, jj2
          if ((i1 == i .and. j1 == j)) cycle  ! Skip self-pair
          if (precep_back(i1,j1) >= 0.0d0 .and. precep_obs(i1,j1) >= 0.0d0) then
            lata = lat1(i1,j1)
            lona = lon1(i1,j1)
            dist = great_circle_distance(latb, lonb, lata, lona)
            diffij1(i,j)=precep_obs(i1,j1)-precep_back(i1,j1)
 !            if (dist <= max_distance_km * 1000.0d0) then
              index = int(dist * 0.001d0 / vario_bin_dist) + 1
              if (index <= max_vario_bins) then
                vario_d(index,nvar_day) = vario_d(index,nvar_day) + &
                (diffij1(i,j)-diffij(i,j))*(diffij1(i,j)-diffij(i,j))
                icounts_vario_d(index,nvar_day) = icounts_vario_d(index,nvar_day) + 1
              end if
!            end if
          end if
        end do
      end do
    end if
  end do
end do

  if (nvar_day>28) then
  write(filenamevar, '(A,I4.4,I2.2,I2.2,A)') 'merra2_semivar/Empvmerr_', actual_year, imonth, iday, '.txt'

  open(unit=50, file=trim(filenamevar), status='replace', action='write', iostat=i)
  if (i /= 0) then
    print *, 'Error opening output file:', trim(filenamevar)
    stop
  end if
 endif

! Aggregate semivariogram over last 28 days
do jj = 1, max_vario_bins
  if (nvar_day <= 28) then
    d1 = 1
  else
    d1 = nvar_day - 28
  end if

  vario(jj) = SUM(vario_d(jj, d1:nvar_day)) 
  icounts_vario(jj) = SUM(icounts_vario_d(jj, d1:nvar_day))

  if (icounts_vario(jj) > 0) then
    vario(jj) = 0.5 * vario(jj) / icounts_vario(jj)
  end if
  if (nvar_day>=28) then
    write(50, '(A,f10.0,A,e13.7,A,I14)') &
      ' dist: ', real((jj-1)*vario_bin_dist), &
      ' semivariogram: ', vario(jj), ' icount: ', icounts_vario(jj)

  end if
end do
close(50)

      end do
    end do
  end do

contains

  subroutine check(status)
    integer, intent(in) :: status
    if (status /= nf90_noerr) then
      print *, trim(nf90_strerror(status))
      stop "Stopped"
    end if
  end subroutine check

subroutine read_obs_data(filename, nlines, lat_obsd, lon_obsd, precepobsd, nobsd1)
    implicit none
    character(len=*), intent(in) :: filename
    integer, intent(out) :: nlines, nobsd1
    real(8), intent(out) :: lat_obsd(:), lon_obsd(:), precepobsd(:)

    integer :: i, unit, ios
    real(8) :: lat, lon, apcp
    character(len=20) :: station_id

    unit = 20
    open(unit, file=filename, status='old', action='read', iostat=ios)
    if (ios /= 0) then
        print *, 'Error opening file: ', trim(filename)
        stop
    end if

    ! First line: number of stations
    read(unit, *, iostat=ios) nlines
    if (ios /= 0) then
        print *, 'Error reading header in: ', trim(filename)
        stop
    end if
    nobsd1 = nlines

    ! Read station data
    do i = 1, nlines
        read(unit, *, iostat=ios) station_id, lat, lon, apcp
        if (ios /= 0) then
            print *, 'Error reading line ', i, ' in file: ', trim(filename)
            stop
        end if
        lat_obsd(i) = lat
        lon_obsd(i) = lon
        precepobsd(i) = apcp
    end do

    close(unit)
end subroutine read_obs_data
   

   !---------------------------------------------------------------------------
   ! Calculates great circle distance between two lat/lon points on globe
   ! using Vincenty formula.
   ! See https://en.wikipedia.org/wiki/Great-circle_distance
   real*8 function great_circle_distance(lat1,lon1,lat2,lon2)

      ! Defaults
      implicit none

      ! Arguments
      real*8, intent(in) :: lat1, lon1, lat2, lon2

      ! Local variables
      real*8 :: radlat1, radlon1, radlat2, radlon2
      real*8 :: pi,deg2rad,lon_abs_diff,central_angle,term1, term2, term3

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

      end
