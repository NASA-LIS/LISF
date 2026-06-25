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
! !ROUTINE: temporal_disaggregation_1km
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
!    temporal_disaggregation_1km.f90 -o temporal_disaggregation_1km 
!    -L/discover/nobackup/projects/lis/libs/sles-12.3/netcdf/4.8.1_intel-2021.4.0/lib 
!    -L/discover/nobackup/projects/lis/libs/sles-12.3/hdf5/1.12.1_intel-2021.4.0/lib 
!    -lnetcdff -lnetcdf -lhdf5_hl -lhdf5 -lz -lcurl 
!    -Wl,--no-relax -shared-intel
!
! !DESCRIPTION:
!  This program temporally disaggregates daily 1-km NLDAS-3 precipitation
!  fields into hourly precipitation using hourly MERRA-2 and IMERG fractions.
!  The output is a daily NetCDF file containing 24 hourly precipitation fields
!  on the 1-km NLDAS-3 grid.
!
!EOP

Program temporal_disaggregation_1km
use netcdf
  implicit none

! ============================================================
! USER SETTINGS
! ============================================================
      integer, parameter :: start_year  = 2024
      integer, parameter :: end_year    = 2024
      integer, parameter :: start_month = 1
      integer, parameter :: end_month   = 2

      character(len=*), parameter :: nldas3_daily_root  = "./nldas3_1km_daily"
      character(len=*), parameter :: merra2_hourly_root = "./merra2_hourly/SURFACEMODEL"
      character(len=*), parameter :: imerg_hourly_root  = "./imerg_hourly/SURFACEMODEL"
      character(len=*), parameter :: output_root        = "./nldas3_1km_hourly"
      character(len=*), parameter :: mask_file          = "lis_input.nldas3.noahmp401.1km.nc"

      character(len=*), parameter :: file_prefix        = "LIS_HIST_"

! ============================================================
! ORIGINAL GRID SETTINGS
! ============================================================
      integer, parameter :: NX=2926, NY=1626, NX1=11700, NY1=6500
      integer, parameter :: ndims=2, ndims1=3, ndt=24

! ============================================================
! ORIGINAL VARIABLE DECLARATIONS
! ============================================================
      integer :: i,j,iyear,imonth,iday,ihour,ihoura
      integer :: nday(12),xi,yi
      integer :: dd

      character(len=300) :: filename1,filename2,filename4
      character(len=300) :: fileplacea,filename1a
      character(len=300) :: daily_filename,hourly_filename
      character(len=4)   :: yyyy
      character(len=2)   :: mm,day_string,hour_string
      character(len=6)   :: yyyymm

      integer :: dimids(ndims),numLons,numLats,dimid,dimids1(ndims1)
      integer :: ncid,totprecip_id,lat_dimid,lon_dimid,lat_id,lon_id
      integer :: latitude_id,longitude_id
      integer :: x,y,ndimsp,nvarsp,nattsp,unlimdimidp,dt_dimid,dt_id

      real*8 :: convergeThresh,inter

      real*8, allocatable, dimension (:,:) :: lat1,lon1
      real*8, allocatable, dimension (:,:) :: precepmerra2,precepimerg,precep_hr
      real*8, allocatable, dimension (:,:) :: mask,precept1,precept2
      real*8, allocatable, dimension (:,:,:) :: preceph1,preceph2,preceph
      real*8, allocatable, dimension (:) :: lat,lon,t_dt

      integer, allocatable, dimension (:) :: start,count

! ============================================================
! ALLOCATE ARRAYS
! ============================================================
      allocate(start(ndims),count(ndims))
      allocate(precepmerra2(NX,NY),precepimerg(NX,NY),precep_hr(NX1,NY1),t_dt(24))
      allocate(precept1(NX,NY),precept2(NX,NY))
      allocate(preceph1(24,NX,NY),preceph2(24,NX,NY),preceph(NX1,NY1,24))
      allocate(lat1(NX1,NY1),lon1(NX1,NY1),mask(NX1,NY1))
      allocate(lat(NY1),lon(NX1))

! ============================================================
! TIME DIMENSION
! ============================================================
      do i=1,24
          t_dt(i)=i
      end do

      preceph=0.0d0

! ============================================================
! DAYS PER MONTH
! ============================================================
      nday(1)=31   ;   nday(2)=28   ;   nday(3)=31   ;   nday(4)=30
      nday(5)=31   ;   nday(6)=30   ;   nday(7)=31   ;   nday(8)=31
      nday(9)=30   ;   nday(10)=31  ;   nday(11)=30  ;   nday(12)=31

! ============================================================
! BUILD 1 KM GRID
! ============================================================
      do i=1,NX1
          do j=1,NY1
              lat1(i,j)=7.005d0+(j-1)*0.01d0
              lat(j)=7.005d0+(j-1)*0.01d0
              lon1(i,j)=-168.995d0+(i-1)*0.01d0
              lon(i)=-168.995d0+(i-1)*0.01d0
          end do
      end do

! ============================================================
! READ 1 KM LANDMASK
! ============================================================
      filename2=trim(mask_file)

      write(*,*) "Reading mask: ", trim(filename2)

      call check(NF90_OPEN(trim(filename2),NF90_NOWRITE,ncid))
      call check(nf90_inquire(ncid,ndimsp,nvarsp,nattsp,unlimdimidp))
      dimid=1
      call check(nf90_inq_varid(ncid,"LANDMASK",totprecip_id))
      call check(nf90_inquire_variable(ncid,totprecip_id,dimids=dimIDs))
      call check(nf90_inquire_dimension(ncid,dimIDs(1),len=numLons))
      call check(nf90_inquire_dimension(ncid,dimIDs(2),len=numLats))
      dimids=(/numLons,numLats/)
      call check(nf90_get_var(ncid,totprecip_id,mask))
      call check(nf90_close(ncid))

      dd=0

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

! ============================================================
! MONTH LOOP
! ============================================================
          DO imonth=start_month,end_month

              WRITE(mm,'(I2.2)') imonth
              yyyymm = yyyy//mm

! ============================================================
! DAY LOOP
! ============================================================
              DO iday=1,nday(imonth)

                  WRITE(day_string,'(I2.2)') iday

                  write(*,*) iyear,imonth,iday

! ============================================================
! READ DAILY 1 KM NLDAS3 PRECIPITATION
! ============================================================
                  daily_filename=trim(file_prefix)//yyyymm//day_string//'0000.d01.nc'

                  fileplacea=trim(nldas3_daily_root)//'/'//yyyymm
                  filename1a=trim(fileplacea)//'/'//trim(daily_filename)

                  write(*,*) "Reading daily NLDAS3 1km: ", trim(filename1a)

                  call check(NF90_OPEN(trim(filename1a), NF90_NOWRITE, ncid))
                  call check(nf90_inquire(ncid, ndimsp, nvarsp, nattsp, unlimdimidp))
                  dimid=1
                  call check(nf90_inq_varid(ncid,"TotalPrecip_tavg",totprecip_id))
                  call check(nf90_inquire_variable(ncid,totprecip_id,dimids = dimIDs))
                  call check(nf90_inquire_dimension(ncid,dimIDs(1),len = numLons))
                  call check(nf90_inquire_dimension(ncid,dimIDs(2),len = numLats))
                  dimids = (/numLons,numLats/)
                  call check(nf90_get_var(ncid,totprecip_id,precep_hr))
                  call check(nf90_close(ncid))

! ============================================================
! CHECK DAILY 1 KM FIELD OVER LAND
! ============================================================
                  do i=1,NX1
                      do j=1,NY1
                          if (mask(i,j)==1.0d0 .and. precep_hr(i,j)<0.0d0) then
                              write(*,*) i,j,precep_hr(i,j)
                          end if
                      end do
                  end do

                  precept1=0.0d0
                  precept2=0.0d0
                  preceph1=0.0d0
                  preceph2=0.0d0

! ============================================================
! HOUR LOOP
! ============================================================
                  DO ihoura=1,24

                      ihour=ihoura-1
                      WRITE(hour_string,'(I2.2)') ihour

                      hourly_filename=trim(file_prefix)//yyyymm//day_string//hour_string//'00.d01.nc'

! ============================================================
! READ HOURLY MERRA2 4 KM
! ============================================================
                      fileplacea=trim(merra2_hourly_root)//'/'//yyyymm
                      filename1a=trim(fileplacea)//'/'//trim(hourly_filename)

                      call check(NF90_OPEN(trim(filename1a), NF90_NOWRITE, ncid))
                      call check(nf90_inquire(ncid, ndimsp, nvarsp, nattsp, unlimdimidp))
                      dimid=1
                      call check(nf90_inq_varid(ncid, "TotalPrecip_tavg", totprecip_id))
                      call check(nf90_inquire_variable(ncid, totprecip_id, dimids = dimIDs))
                      call check(nf90_inquire_dimension(ncid, dimIDs(1), len = numLons))
                      call check(nf90_inquire_dimension(ncid, dimIDs(2), len = numLats))
                      dimids = (/numLons, numLats /)
                      call check(nf90_get_var(ncid, totprecip_id, precepmerra2))
                      call check(nf90_close(ncid))

! ============================================================
! READ HOURLY IMERG 4 KM
! ============================================================
                      fileplacea=trim(imerg_hourly_root)//'/'//yyyymm
                      filename2=trim(fileplacea)//'/'//trim(hourly_filename)

                      call check(NF90_OPEN(trim(filename2), NF90_NOWRITE, ncid))
                      call check(nf90_inquire(ncid, ndimsp, nvarsp, nattsp, unlimdimidp))
                      dimid=1
                      call check(nf90_inq_varid(ncid, "TotalPrecip_tavg", totprecip_id))
                      call check(nf90_inquire_variable(ncid, totprecip_id, dimids = dimIDs))
                      call check(nf90_inquire_dimension(ncid, dimIDs(1), len = numLons))
                      call check(nf90_inquire_dimension(ncid, dimIDs(2), len = numLats))
                      dimids = (/numLons, numLats /)
                      call check(nf90_get_var(ncid, totprecip_id, precepimerg))
                      call check(nf90_close(ncid))

! ============================================================
! ACCUMULATE DAILY TOTALS AND STORE HOURLY VALUES
! ============================================================
                      do i=1,NX
                          do j=1,NY
                              precept1(i,j)=precept1(i,j)+precepmerra2(i,j)
                              preceph1(ihoura,i,j)=precepmerra2(i,j)

                              precept2(i,j)=precept2(i,j)+precepimerg(i,j)
                              preceph2(ihoura,i,j)=precepimerg(i,j)
                          end do
                      end do

                  END DO   ! HOUR

! ============================================================
! TEMPORAL DISAGGREGATION TO 24 HOURS AT 1 KM
! ============================================================
                  preceph=0.0d0

                  DO ihoura=1,24

                      inter=0.0d0

                      do i=1,NX1
                          do j=1,NY1

                              xi=int((lon1(i,j)-(-169.0d0))/0.04d0)+1
                              yi=int((lat1(i,j)-7.0d0)/0.04d0)+1

                              IF (precept1(xi,yi)>0.0d0) THEN
                                  inter=preceph1(ihoura,xi,yi)/precept1(xi,yi)
                              ELSE IF (precept2(xi,yi)>0.0d0) THEN
                                  inter=preceph2(ihoura,xi,yi)/precept2(xi,yi)
                              ELSE
                                  inter=1.0d0/24.0d0
                              END IF

                              if(precep_hr(i,j)>=0.0d0) then
                                  preceph(i,j,ihoura)=precep_hr(i,j)*inter
                              else
                                  preceph(i,j,ihoura)=-9999.0d0
                              end if

                              if (mask(i,j)==1.0d0 .and. preceph(i,j,ihoura)==-9999.0d0) then
                                  if (precep_hr(i,j) .ne. -9999.0d0) then
                                      preceph(i,j,ihoura)=precep_hr(i,j)/24.0d0
                                  else
                                      write(*,*) 'notfixed',i,j
                                  end if
                              end if

                              if (mask(i,j)==1.0d0 .and. preceph(i,j,ihoura)<0.0d0) then
                                  write(*,*) i,j,preceph(i,j,ihoura),precep_hr(i,j)
                              end if

                          end do
                      end do

                  END DO   ! HOUR

! ============================================================
! WRITE DAILY FILE WITH 24 HOURLY FIELDS
! ============================================================
                  filename1=trim(daily_filename)

                  fileplacea=trim(output_root)//'/'//yyyymm
                  call execute_command_line("mkdir -p " // trim(fileplacea))

                  filename1a=trim(fileplacea)//'/'//trim(filename1)

                  write(*,*) "Writing hourly 1km output: ", trim(filename1a)

                  call check(NF90_create(trim(filename1a),NF90_clobber,ncid))

                  call check(nf90_def_dim(ncid,"latitude",NY1,lat_dimid))
                  call check(nf90_def_dim(ncid,"longitude",NX1,lon_dimid))

                  call check(nf90_def_var(ncid,"latitude",NF90_REAL,lat_dimid,lat_id))
                  call check(nf90_def_var(ncid,"longitude",NF90_REAL,lon_dimid,lon_id))

                  dimids=(/lon_dimid,lat_dimid/)

                  call check(nf90_def_var(ncid,"lat",NF90_REAL,dimids,latitude_id))
                  CALL check(nf90_put_att(ncid,latitude_id,"units","degree_north"))
                  call check(nf90_put_att(ncid,latitude_id,"_FillValue", -9999.00))
                  CALL check(nf90_put_att(ncid,latitude_id,"long_name","latitude"))

                  call check(nf90_def_var(ncid,"lon",NF90_REAL,dimids,longitude_id))
                  CALL check(nf90_put_att(ncid,longitude_id,"units","degree_east"))
                  call check(nf90_put_att(ncid,longitude_id,"_FillValue", -9999.00))
                  CALL check(nf90_put_att(ncid,longitude_id,"long_name","longitude"))

                  call check(nf90_def_dim(ncid, "TotalPrecip_hours", Ndt, dt_dimid))
                  call check(nf90_def_var(ncid, "TotalPrecip_hours", NF90_REAL, dt_dimid, dt_id))

                  dimids1=(/lon_dimid,lat_dimid,dt_dimid/)

                  call check(nf90_def_var(ncid,"TotalPrecip",NF90_REAL,dimids1,totprecip_id))
                  CALL check(nf90_put_att(ncid,totprecip_id,"units","mm"))
                  call check(nf90_put_att(ncid,totprecip_id,"_FillValue", -9999.00))
                  CALL check(nf90_put_att(ncid,totprecip_id,"long_name","total precipitation"))

                  CALL check(nf90_enddef(ncid))

                  call check(nf90_put_var(ncid,lat_id,lat))
                  call check(nf90_put_var(ncid,lon_id,lon))
                  call check(nf90_put_var(ncid,dt_id,t_dt))
                  call check(nf90_put_var(ncid,latitude_id,lat1))
                  call check(nf90_put_var(ncid,longitude_id,lon1))
                  call check(nf90_put_var(ncid,totprecip_id,preceph))

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

end Program temporal_disaggregation_1km
