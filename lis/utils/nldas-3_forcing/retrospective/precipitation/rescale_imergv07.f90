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
! !ROUTINE: rescale_imergv07
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
!    rescale_imergv07.f90 -o rescale_imergv07 
!    -L/discover/nobackup/projects/lis/libs/sles-12.3/netcdf/4.8.1_intel-2021.4.0/lib 
!    -L/discover/nobackup/projects/lis/libs/sles-12.3/hdf5/1.12.1_intel-2021.4.0/lib 
!    -lnetcdff -lnetcdf -lhdf5_hl -lhdf5 -lz -lcurl 
!    -Wl,--no-relax -shared-intel
!
! !DESCRIPTION:
!  This program rescales daily IMERG V07 precipitation fields using monthly
!  scale factors from GPCP. The output is a daily NetCDF precipitation file on the
!  NLDAS-3 grid.
!
!EOP

Program rescale_imergv07
use netcdf
  implicit none

! ============================================================
! USER SETTINGS
! Change only this block
! ============================================================
      integer, parameter :: start_year  = 2026
      integer, parameter :: end_year    = 2026
      integer, parameter :: start_month = 1
      integer, parameter :: end_month   = 12

      character(len=*), parameter :: scale_input_dir    = "./gimerg_month"
      character(len=*), parameter :: imerg_input_root   = "./imerg10km_v07/SURFACEMODEL"
      character(len=*), parameter :: output_root        = "./gimerg10km_daily"

      character(len=*), parameter :: scale_prefix       = "GPERM_CLIM_"
      character(len=*), parameter :: daily_prefix       = "LIS_HIST_"

! ============================================================
! ORIGINAL GRID SETTINGS
! ============================================================
      integer, parameter :: NX=1171, NY=651, NX1=234, NY1=130, NDT=1

! ============================================================
! ORIGINAL VARIABLE DECLARATIONS
! ============================================================
      integer :: i,j,iyear,imonth,iday
      integer :: nday(12)
      integer :: x,y
      integer :: inv_p
      integer :: dday

      Character (LEN=300) :: filename1,filename2,filename4
      Character (LEN=300) :: fileplace1,fileplacea,filename1a
      Character (LEN=4)   :: yyyy
      Character (LEN=2)   :: yy,mm,dd
      Character (LEN=6)   :: yyyymm

      Integer, parameter :: ndims=2,ndims1=2
      integer :: dimids(ndims),numLons,numLats,dimid,dimids1(ndims1),numndts

      real*8 :: convergeThresh,lata,lona,latb,lonb,dx,dx1

      real, allocatable, dimension (:,:) :: latitude,longitude,latitude1,longitude1
      real*8, allocatable, dimension (:,:) :: gimerg,scal,imerg_v07
      real*8, allocatable, dimension (:,:,:) :: gpcp
      real*8, allocatable, dimension (:) :: lat,lon

      integer :: ncid, totprecip_id,lat_dimid,lon_dimid,lat_id,lon_id
      integer :: latitude_id,longitude_id
      integer :: ndimsp,nvarsp,nattsp,unlimdimidp,cfm_id

      integer, allocatable, dimension (:) :: start,count
      integer, allocatable, dimension (:,:) :: id_lr,jd_lr,cell,cell1

! ============================================================
! ALLOCATE ARRAYS
! ============================================================
      allocate(start(ndims),count(ndims))
      allocate(gimerg(NX,NY),scal(NX,NY),imerg_v07(NX,NY),gpcp(NX1,NY1,1))
      allocate(latitude(NX,NY),longitude(NX,NY))
      allocate(lat(NY),lon(NX))
      allocate(latitude1(NX1,NY1),longitude1(NX1,NY1))

! ============================================================
! DAYS PER MONTH
! ============================================================
      nday(1)=31   ;   nday(2)=28   ;   nday(3)=31   ;   nday(4)=30
      nday(5)=31   ;   nday(6)=30   ;   nday(7)=31   ;   nday(8)=31
      nday(9)=30   ;   nday(10)=31  ;   nday(11)=30  ;   nday(12)=31

      inv_p=2
      dx=0.1d0
      dx1=0.5d0

! ============================================================
! BUILD 0.5 DEGREE GRID
! ============================================================
      do i=1,NX1
          do j=1,NY1
              latitude1(i,j)=7.0+(dx1*(j-1))
              longitude1(i,j)=-169+(dx1*(i-1))
          end do
      end do

! ============================================================
! BUILD 0.1 DEGREE GRID
! ============================================================
      do i=1,NX
          do j=1,NY
              latitude(i,j)=7.0+(dx*(j-1))
              longitude(i,j)=-169+(dx*(i-1))
              lat(j)=7.0+(dx*(j-1))
              lon(i)=-169.0+(dx*(i-1))
          end do
      end do

! ============================================================
! YEAR LOOP
! ============================================================
      DO iyear=start_year,end_year

          IF (is_leap_year(iyear)) THEN
              nday(2)=29
          ELSE
              nday(2)=28
          END IF

          WRITE(yyyy,'(I4.4)') iyear
          yy = yyyy(3:4)

! ============================================================
! MONTH LOOP
! ============================================================
          DO imonth=start_month,end_month

              WRITE(mm,'(I2.2)') imonth
              yyyymm = yyyy//mm

! ============================================================
! READ MONTHLY SCALE FACTOR
! ============================================================
              filename2=trim(scale_input_dir)//'/'//trim(scale_prefix)//mm//'.nc'

              write(*,*) "Reading scale file: ", trim(filename2)

              call check(NF90_OPEN(trim(filename2), NF90_NOWRITE, ncid))
              call check(nf90_inquire(ncid, ndimsp, nvarsp, nattsp, unlimdimidp))
              dimid=1
              call check(nf90_inq_varid(ncid, "Scale", totprecip_id))
              call check(nf90_inquire_variable(ncid, totprecip_id, dimids = dimIDs))
              call check(nf90_inquire_dimension(ncid, dimIDs(1), len = numLons))
              call check(nf90_inquire_dimension(ncid, dimIDs(2), len = numLats))
              dimids = (/numLons, numLats /)
              call check(nf90_get_var(ncid, totprecip_id,scal))
              call check(nf90_close(ncid))

              write(*,*) iyear,imonth

! ============================================================
! DAY LOOP
! ============================================================
              DO iday=1,nday(imonth)

                  WRITE(dd,'(I2.2)') iday

! ============================================================
! DAILY IMERG FILE NAME
! ============================================================
                  filename1=trim(daily_prefix)//yyyymm//dd//'0000.d01.nc'

                  fileplacea=trim(imerg_input_root)//'/'//yyyymm
                  filename1a=trim(fileplacea)//'/'//trim(filename1)

                  write(*,*) "Reading IMERG: ", trim(filename1a)

                  call check(NF90_OPEN(trim(filename1a), NF90_NOWRITE, ncid))
                  call check(nf90_inquire(ncid, ndimsp, nvarsp, nattsp, unlimdimidp))
                  dimid=1
                  call check(nf90_inq_varid(ncid, "TotalPrecip_tavg", totprecip_id))
                  call check(nf90_inquire_variable(ncid, totprecip_id, dimids = dimIDs))
                  call check(nf90_inquire_dimension(ncid, dimIDs(1), len = numLons))
                  call check(nf90_inquire_dimension(ncid, dimIDs(2), len = numLats))
                  dimids = (/numLons, numLats /)
                  call check(nf90_get_var(ncid, totprecip_id,imerg_v07))
                  call check(nf90_close(ncid))

! ============================================================
! APPLY MONTHLY SCALE FACTOR
! ============================================================
                  do i=1,NX
                      do j=1,NY
                          if (scal(i,j)<0.0d0) then
                              gimerg(i,j)=imerg_v07(i,j)
                          else
                              gimerg(i,j)=imerg_v07(i,j)*scal(i,j)
                          end if
                      end do
                  end do

! ============================================================
! WRITE DAILY OUTPUT
! ============================================================
                  fileplace1=trim(output_root)//'/'//yyyymm
                  call execute_command_line("mkdir -p " // trim(fileplace1))

                  filename4=trim(fileplace1)//'/'//trim(filename1)

                  write(*,*) "Writing output: ", trim(filename4)

                  call check(NF90_create(trim(filename4),NF90_clobber,ncid))

                  call check(nf90_def_dim(ncid,"latitude",NY,lat_dimid))
                  call check(nf90_def_dim(ncid,"longitude",NX,lon_dimid))

                  call check(nf90_def_var(ncid,"latitude",NF90_REAL,lat_dimid,lat_id))
                  call check(nf90_def_var(ncid,"longitude",NF90_REAL,lon_dimid,lon_id))

                  dimids=(/lon_dimid,lat_dimid/)

                  call check(nf90_def_var(ncid,"TotalPrecip",NF90_REAL,dimids,totprecip_id))
                  call check(nf90_def_var(ncid,"lat",NF90_REAL,dimids,latitude_id))
                  call check(nf90_def_var(ncid,"lon",NF90_REAL,dimids,longitude_id))

                  CALL check(nf90_enddef(ncid))

                  call check(nf90_put_var(ncid,lat_id,lat))
                  call check(nf90_put_var(ncid,lon_id,lon))
                  call check(nf90_put_var(ncid,totprecip_id,gimerg))
                  call check(nf90_put_var(ncid,latitude_id,latitude))
                  call check(nf90_put_var(ncid,longitude_id,longitude))

                  CALL check(nf90_close(ncid))

              END DO   ! DAY

          END DO   ! MONTH

      END DO   ! YEAR

 14 FORMAT(2000(1x,e13.7))

contains

  subroutine check(status)
    integer, intent(in) :: status

    if(status /= nf90_noerr) then
      print *, trim(nf90_strerror(status))
      stop "Stopped"
    end if
  end subroutine check


  logical function is_leap_year(year)
    implicit none
    integer, intent(in) :: year

    is_leap_year = .false.

    if (mod(year,400) == 0) then
        is_leap_year = .true.
    else if (mod(year,100) == 0) then
        is_leap_year = .false.
    else if (mod(year,4) == 0) then
        is_leap_year = .true.
    end if

  end function is_leap_year


   !---------------------------------------------------------------------------
   ! Kept from original code, although it is not used in this program.
   !---------------------------------------------------------------------------
   real*8 function great_circle_distance(lat1,lon1,lat2,lon2)

      implicit none

      real*8, intent(in) :: lat1, lon1, lat2, lon2

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

end Program rescale_imergv07
