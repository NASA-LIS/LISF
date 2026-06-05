Program downscaled_merra2_cloud_4km
use netcdf
  implicit none

! ============================================================
! USER SETTINGS
! Change only this block
! ============================================================
      integer, parameter :: start_year  = 2026
      integer, parameter :: end_year    = 2026
      integer, parameter :: start_month = 3
      integer, parameter :: end_month   = 12

      character(len=*), parameter :: modis_cfm_dir       = "./Cloudcover_4km"
      character(len=*), parameter :: myd_cfm_dir         = "./Cloudcover_myd_4km"
      character(len=*), parameter :: merra2_input_root   = "./merra2_daily/SURFACEMODEL"
      character(len=*), parameter :: merra2_output_root  = "./merra2_daily_4km"

! ============================================================
! ORIGINAL GRID SETTINGS
! ============================================================
      integer, parameter :: NX=2926, NY=1626, NX1=226, NY1=126

! ============================================================
! ORIGINAL VARIABLE DECLARATIONS
! ============================================================
      integer :: i,j,iyear,imonth,iday
      integer :: nday(12),ii1,ii2,jj1,jj2,xi,yi
      integer :: x,y
      integer :: inv_p

      Character (LEN=300) :: filename1,filename2,filename4
      Character (LEN=300) :: fileplacea,fileplace1,filename1a
      Character (LEN=4)   :: yyyy
      Character (LEN=2)   :: mm,dd
      Character (LEN=6)   :: yyyymm

      Integer, parameter :: ndims=2, ndims1=2
      integer :: dimids(ndims),numLons,numLats,dimid,dimids1(ndims1)

      real*8 :: dist,inv_dist,inv_dist_pt,merra2_pt,cfm_pt
      real*8 :: lata,lona,latb,lonb,dx,dx1,merra2_max

      real, allocatable, dimension (:,:) :: latitude,longitude,latitude1,longitude1

      real*8, allocatable, dimension (:,:) :: merra2_hr
      real*8, allocatable, dimension (:,:) :: cfm_50km,inv_dist_50km,pcfm_50km
      real*8, allocatable, dimension (:,:) :: merra2_4km,cfm_4km,merra2_50km
      real*8, allocatable, dimension (:,:) :: cfm_4km_mod,cfm_4km_myd

      real*8, allocatable, dimension (:) :: lat,lon

      integer :: ncid, totprecip_id,lat_dimid,lon_dimid,lat_id,lon_id
      integer :: latitude_id,longitude_id
      integer :: ndimsp,nvarsp,nattsp,unlimdimidp,cfm_id

! ============================================================
! ALLOCATE ARRAYS
! ============================================================
      allocate(cfm_50km(NX1,NY1),cfm_4km(NX,NY),inv_dist_50km(NX1,NY1))
      allocate(latitude(NX,NY),longitude(NX,NY),merra2_hr(NX,NY))
      allocate(lat(NY),lon(NX),merra2_4km(NX,NY))
      allocate(cfm_4km_mod(NX,NY),cfm_4km_myd(NX,NY),pcfm_50km(NX1,NY1))
      allocate(latitude1(NX1,NY1),longitude1(NX1,NY1),merra2_50km(NX1,NY1))

! ============================================================
! DAYS PER MONTH
! February is updated inside the year loop
! ============================================================
      nday(1)=31   ;   nday(2)=28   ;   nday(3)=31   ;   nday(4)=30
      nday(5)=31   ;   nday(6)=30   ;   nday(7)=31   ;   nday(8)=31
      nday(9)=30   ;   nday(10)=31  ;   nday(11)=30  ;   nday(12)=31

      inv_p=2
      dx=0.04d0
      dx1=0.52d0

! ============================================================
! BUILD 50 KM GRID
! ============================================================
      do i=1,NX1
          do j=1,NY1
              latitude1(i,j)=7.0+(dx1*(j-1))
              longitude1(i,j)=-169+(dx1*(i-1))
          end do
      end do

! ============================================================
! BUILD 4 KM GRID
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
! Same structure as original, but year range is user-defined
! ============================================================
      DO iyear=start_year,end_year

! ============================================================
! LEAP YEAR HANDLING
! ============================================================
          IF (is_leap_year(iyear)) THEN
              nday(2)=29
          ELSE
              nday(2)=28
          END IF

          WRITE(yyyy,'(I4.4)') iyear

! ============================================================
! MONTH LOOP
! Same structure as original, but month range is user-defined
! ============================================================
          DO imonth=start_month,end_month

              WRITE(mm,'(I2.2)') imonth
              yyyymm = yyyy//mm

! ============================================================
! READ MONTHLY MODIS CFM 4 KM
! ============================================================
              filename1='CFM_'//yyyymm//'.nc'
              filename2=trim(modis_cfm_dir)//'/'//trim(filename1)

              write(*,*) "Reading MODIS CFM: ", trim(filename2)

              call check(NF90_OPEN(trim(filename2),NF90_NOWRITE,ncid))
              call check(nf90_inquire(ncid,ndimsp,nvarsp,nattsp,unlimdimidp))
              dimid=1
              call check(nf90_inq_varid(ncid,"CFM",totprecip_id))
              call check(nf90_inquire_variable(ncid,totprecip_id,dimids=dimIDs))
              call check(nf90_inquire_dimension(ncid,dimIDs(1),len=numLons))
              call check(nf90_inquire_dimension(ncid,dimIDs(2),len=numLats))
              dimids=(/numLons,numLats/)
              call check(nf90_get_var(ncid,totprecip_id,cfm_4km_mod))
              call check(nf90_close(ncid))

! ============================================================
! READ MONTHLY MYD CFM 4 KM
! ============================================================
              filename2=trim(myd_cfm_dir)//'/'//trim(filename1)

              write(*,*) "Reading MYD CFM: ", trim(filename2)

              call check(NF90_OPEN(trim(filename2),NF90_NOWRITE,ncid))
              call check(nf90_inquire(ncid,ndimsp,nvarsp,nattsp,unlimdimidp))
              dimid=1
              call check(nf90_inq_varid(ncid,"CFM",totprecip_id))
              call check(nf90_inquire_variable(ncid,totprecip_id,dimids=dimIDs))
              call check(nf90_inquire_dimension(ncid,dimIDs(1),len=numLons))
              call check(nf90_inquire_dimension(ncid,dimIDs(2),len=numLats))
              dimids=(/numLons,numLats/)
              call check(nf90_get_var(ncid,totprecip_id,cfm_4km_myd))
              call check(nf90_close(ncid))

! ============================================================
! AVERAGE MODIS AND MYD CFM
! ============================================================
              do x=1,NX
                  do y=1,NY
                      cfm_4km(x,y)=(cfm_4km_myd(x,y)+cfm_4km_mod(x,y))/2.0d0
                  end do
              end do

! ============================================================
! AGGREGATE CFM FROM 4 KM TO 50 KM
! ============================================================
              cfm_50km=0.0d0
              inv_dist_50km=0.0d0

              do i=1,NX1
                  do j=1,NY1

                      xi=int((longitude1(i,j)-(-169.0d0))/dx)+1
                      yi=int((latitude1(i,j)-7.0d0)/dx)+1

                      lata=latitude1(i,j)
                      lona=longitude1(i,j)

                      ii1=xi-13
                      ii2=xi+13
                      jj1=yi-13
                      jj2=yi+13

                      if (ii1<1) ii1=1
                      if (ii1>NX) ii1=NX
                      if (ii2>NX) ii2=NX
                      if (ii2<1) ii2=1

                      if (jj1<1) jj1=1
                      if (jj1>NY) jj1=NY
                      if (jj2>NY) jj2=NY
                      if (jj2<1) jj2=1

                      do x=ii1,ii2
                          do y=jj1,jj2

                              if (cfm_4km(x,y)>=0.0d0) then
                                  if (abs(latitude1(i,j)-latitude(x,y))<=dx1) then
                                      if (abs(abs(longitude1(i,j))-abs(longitude(x,y)))<=dx1) then

                                          latb=latitude(x,y)
                                          lonb=longitude(x,y)

                                          dist=great_circle_distance(latb,lonb,lata,lona)

                                          if (dist==0.0d0) then
                                              inv_dist=1.0d0
                                          else
                                              dist=(dist/1000.0d0)*0.01d0
                                              inv_dist=(1.0d0/dist)**inv_p
                                          end if

                                          inv_dist_50km(i,j)=inv_dist_50km(i,j)+inv_dist
                                          cfm_50km(i,j)=cfm_50km(i,j)+inv_dist*cfm_4km(x,y)

                                      end if
                                  end if
                              end if

                          end do
                      end do

                      if (inv_dist_50km(i,j)==0.0d0) then
                          cfm_50km(i,j)=0.0d0
                      else
                          cfm_50km(i,j)=cfm_50km(i,j)/inv_dist_50km(i,j)
                      end if

                  end do
              end do

! ============================================================
! DAY LOOP
! ============================================================
              DO iday=1,nday(imonth)

                  WRITE(dd,'(I2.2)') iday

                  write(*,*) iyear,imonth,iday

! ============================================================
! DAILY INPUT FILE
! ============================================================
                  filename1='LIS_HIST_'//yyyymm//dd//'0000.d01.nc'

                  fileplacea=trim(merra2_input_root)//'/'//yyyymm
                  filename1a=trim(fileplacea)//'/'//trim(filename1)

                  write(*,*) "Reading precipitation: ", trim(filename1a)

                  call check(NF90_OPEN(trim(filename1a), NF90_NOWRITE, ncid))
                  call check(nf90_inquire(ncid, ndimsp, nvarsp, nattsp, unlimdimidp))
                  dimid=1
                  call check(nf90_inq_varid(ncid, "TotalPrecip_tavg", totprecip_id))
                  call check(nf90_inquire_variable(ncid, totprecip_id, dimids = dimIDs))
                  call check(nf90_inquire_dimension(ncid, dimIDs(1), len = numLons))
                  call check(nf90_inquire_dimension(ncid, dimIDs(2), len = numLats))
                  dimids = (/numLons, numLats /)
                  call check(nf90_get_var(ncid, totprecip_id,merra2_4km))
                  call check(nf90_close(ncid))

                  merra2_max=MAXVAL(merra2_4km(:,:))

! ============================================================
! AGGREGATE DAILY MERRA2/LIS PRECIPITATION TO 50 KM
! ============================================================
                  merra2_50km=0.0d0
                  inv_dist_50km=0.0d0
                  pcfm_50km=0.0d0

                  do i=1,NX1
                      do j=1,NY1

                          xi=int((longitude1(i,j)-(-169.0d0))/0.04d0)+1
                          yi=int((latitude1(i,j)-7.0d0)/0.04d0)+1

                          lata=latitude1(i,j)
                          lona=longitude1(i,j)

                          ii1=xi-13
                          ii2=xi+13
                          jj1=yi-13
                          jj2=yi+13

                          if (ii1<1) ii1=1
                          if (ii1>NX) ii1=NX
                          if (ii2>NX) ii2=NX
                          if (ii2<1) ii2=1

                          if (jj1<1) jj1=1
                          if (jj1>NY) jj1=NY
                          if (jj2>NY) jj2=NY
                          if (jj2<1) jj2=1

                          do x=ii1,ii2
                              do y=jj1,jj2

                                  if (cfm_4km(x,y)>=0.0d0) then
                                      if (abs(latitude1(i,j)-latitude(x,y))<=dx1) then
                                          if (abs(abs(longitude1(i,j))-abs(longitude(x,y)))<=dx1) then

                                              latb=latitude(x,y)
                                              lonb=longitude(x,y)

                                              dist=great_circle_distance(latb,lonb,lata,lona)

                                              if (dist==0.0d0) then
                                                  inv_dist=1.0d0
                                              else
                                                  dist=(dist/1000.0d0)*0.01d0
                                                  inv_dist=(1.0d0/dist)**inv_p
                                              end if

                                              inv_dist_50km(i,j)=inv_dist_50km(i,j)+inv_dist
                                              merra2_50km(i,j)=merra2_50km(i,j)+inv_dist*merra2_4km(x,y)
                                              pcfm_50km(i,j)=pcfm_50km(i,j)+inv_dist*merra2_4km(x,y)*cfm_4km(x,y)

                                          end if
                                      end if
                                  end if

                              end do
                          end do

                          if (inv_dist_50km(i,j)==0.0d0) then
                              merra2_50km(i,j)=0.0d0
                              pcfm_50km(i,j)=0.0d0
                          else
                              merra2_50km(i,j)=merra2_50km(i,j)/inv_dist_50km(i,j)
                              pcfm_50km(i,j)=pcfm_50km(i,j)/inv_dist_50km(i,j)
                          end if

                      end do
                  end do

! ============================================================
! APPLY CLOUD-FREQUENCY-BASED DOWNSCALING
! ============================================================
                  do i=1,NX
                      do j=1,NY

                          if (cfm_4km(i,j)>=0.0d0) then

                              inv_dist_pt=0.0d0
                              merra2_pt=0.0d0
                              cfm_pt=0.0d0

                              xi=int((longitude(i,j)-(-169.0d0))/dx1)+1
                              yi=int((latitude(i,j)-7.0d0)/dx1)+1

                              lata=latitude(i,j)
                              lona=longitude(i,j)

                              ii1=xi-1
                              ii2=xi+1
                              jj1=yi-1
                              jj2=yi+1

                              if (ii1<1) ii1=1
                              if (ii1>NX1) ii1=NX1
                              if (ii2>NX1) ii2=NX1
                              if (ii2<1) ii2=1

                              if (jj1<1) jj1=1
                              if (jj1>NY1) jj1=NY1
                              if (jj2>NY1) jj2=NY1
                              if (jj2<1) jj2=1

                              do x=ii1,ii2
                                  do y=jj1,jj2

                                      latb=latitude1(x,y)
                                      lonb=longitude1(x,y)

                                      if (abs(latitude(i,j)-latitude1(x,y))<=dx1) then
                                          if (abs(abs(longitude(i,j))-abs(longitude1(x,y)))<=dx1) then

                                              dist=great_circle_distance(latb,lonb,lata,lona)

                                              if (dist==0.0d0) then
                                                  inv_dist=1.0d0
                                              else
                                                  dist=(dist/1000.0d0)*0.01d0
                                                  inv_dist=(1.0d0/dist)**inv_p
                                              end if

                                              inv_dist_pt=inv_dist_pt+inv_dist
                                              merra2_pt=merra2_pt+inv_dist*merra2_50km(x,y)
                                              cfm_pt=cfm_pt+inv_dist*pcfm_50km(x,y)

                                          end if
                                      end if

                                  end do
                              end do

                              if (inv_dist_pt==0.0d0) then
                                  merra2_pt=0.0d0
                                  cfm_pt=0.0d0
                              else
                                  merra2_pt=merra2_pt/inv_dist_pt
                                  cfm_pt=cfm_pt/inv_dist_pt
                              end if

                              if (cfm_pt==0.0d0) then
                                  merra2_hr(i,j)=merra2_4km(i,j)
                              else
                                  merra2_hr(i,j)=((cfm_4km(i,j)*merra2_4km(i,j))/cfm_pt)*merra2_pt
                              end if

                              if (merra2_hr(i,j)>(1.5d0*merra2_pt)) merra2_hr(i,j)=merra2_4km(i,j)

                          else

                              merra2_hr(i,j)=merra2_4km(i,j)

                          end if

                      end do
                  end do

! ============================================================
! WRITE DAILY OUTPUT
! ============================================================
                  fileplace1=trim(merra2_output_root)//'/'//yyyymm
                  call execute_command_line("mkdir -p " // trim(fileplace1))

                  filename4=trim(fileplace1)//'/'//trim(filename1)

                  write(*,*) "Writing output: ", trim(filename4)

                  call check(NF90_create(trim(filename4),NF90_clobber,ncid))

                  call check(nf90_def_dim(ncid,"latitude",NY,lat_dimid))
                  call check(nf90_def_dim(ncid,"longitude",NX,lon_dimid))

                  call check(nf90_def_var(ncid,"latitude",NF90_REAL,lat_dimid,lat_id))
                  call check(nf90_def_var(ncid,"longitude",NF90_REAL,lon_dimid,lon_id))

                  dimids1=(/lon_dimid,lat_dimid/)

                  call check(nf90_def_var(ncid,"TotalPrecip",NF90_REAL,dimids1,totprecip_id))
!                 call check(nf90_def_var(ncid,"CFM",NF90_REAL,dimids1,cfm_id))

                  call check(nf90_def_var(ncid,"lat",NF90_REAL,dimids1,latitude_id))
                  call check(nf90_def_var(ncid,"lon",NF90_REAL,dimids1,longitude_id))

                  CALL check(nf90_enddef(ncid))

                  call check(nf90_put_var(ncid,lat_id,lat))
                  call check(nf90_put_var(ncid,lon_id,lon))
                  call check(nf90_put_var(ncid,totprecip_id,merra2_hr))
!                 call check(nf90_put_var(ncid,cfm_id,cfm_4km))
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
   ! Calculates great circle distance between two lat/lon points on globe
   ! using Vincenty formula.
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

end Program downscaled_merra2_cloud_4km
