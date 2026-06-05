Program downscaled_nldas3_cloud_1km
use netcdf
  implicit none

! ============================================================
! USER SETTINGS
! Change only this block
! ============================================================
      integer, parameter :: start_year  = 2026
      integer, parameter :: end_year    = 2026
      integer, parameter :: start_month = 1
      integer, parameter :: end_month   = 2

      character(len=*), parameter :: cfm_modis_dir       = "./Cloudcover_1km_real"
      character(len=*), parameter :: cfm_myd_dir         = "./Cloudcover_myd_1km_real"
      character(len=*), parameter :: nldas3_input_root   = "./nldas3_10km_braseth"
      character(len=*), parameter :: merra2_input_root   = "./merra2_daily_4km"
      character(len=*), parameter :: output_root         = "./nldas3_1km_daily"
      character(len=*), parameter :: mask_file           = "lis_input.nldas3.noahmp401.1km.nc"

! ============================================================
! ORIGINAL GRID SETTINGS
! ============================================================
      integer, parameter :: NX=11700, NY=6500, NX1=2926, NY1=1626

! ============================================================
! ORIGINAL VARIABLE DECLARATIONS
! ============================================================
      integer :: i,j,iyear,imonth,iday
      integer :: nday(12),ii1,ii2,jj1,jj2,xi,yi
      integer :: x,y
      integer :: inv_p
      integer :: numndts

      Character (LEN=300) :: filename1,filename2,filename4
      Character (LEN=300) :: fileplacea,fileplace1,filename1a
      Character (LEN=4)   :: yyyy
      Character (LEN=2)   :: mm,dd
      Character (LEN=6)   :: yyyymm

      Integer, parameter :: ndims=2,ndims1=2
      integer :: dimids(ndims),numLons,numLats,dimid,dimids1(ndims1)

      real*8 :: dist,inv_dist,inv_dist_pt,nldas3_pt,cfm_pt
      real*8 :: lata,lona,latb,lonb,dx,dx1,nldas3_hr,inv_dist_4km
      real*8 :: convergeThresh

      real, allocatable, dimension (:,:) :: latitude,longitude,latitude1,longitude1

      real*8, allocatable, dimension (:,:) :: cfm_1km,pcfm_4km,nldas3,mask,merra2_4km
      real*8, allocatable, dimension (:,:) :: nldas3_4km,cfm_4km,nldas3_1km
      real*8, allocatable, dimension (:,:) :: cfm_1km_mod,cfm_1km_myd

      real*8, allocatable, dimension (:) :: lat,lon,t_dt

      integer :: ncid, totprecip_id,lat_dimid,lon_dimid,lat_id,lon_id
      integer :: latitude_id,longitude_id
      integer :: ndimsp,nvarsp,nattsp,unlimdimidp,cfm_id,dt_dimid,dt_id

      integer, allocatable, dimension (:) :: start,count

! ============================================================
! ALLOCATE ARRAYS
! ============================================================
      allocate(start(ndims),count(ndims),cfm_1km(NX,NY))

      allocate(latitude(NX,NY),longitude(NX,NY),nldas3(NX,NY),mask(NX,NY))
      allocate(merra2_4km(NX1,NY1))

      allocate(lat(NY),lon(NX),nldas3_4km(NX1,NY1))
      allocate(cfm_1km_mod(NX,NY),cfm_1km_myd(NX,NY),pcfm_4km(NX1,NY1))

      allocate(latitude1(NX1,NY1),longitude1(NX1,NY1),nldas3_1km(NX,NY))

! ============================================================
! DAYS PER MONTH
! February is updated inside the year loop
! ============================================================
      nday(1)=31   ;   nday(2)=28   ;   nday(3)=31   ;   nday(4)=30
      nday(5)=31   ;   nday(6)=30   ;   nday(7)=31   ;   nday(8)=31
      nday(9)=30   ;   nday(10)=31  ;   nday(11)=30  ;   nday(12)=31

      inv_p=2
      dx=0.01d0
      dx1=0.04d0

! ============================================================
! BUILD 1 KM GRID
! Same as original
! ============================================================
      do i=1,NX
          do j=1,NY
              latitude(i,j)=7.005+(0.01*(j-1))
              lat(j)=7.005+(0.01*(j-1))
              longitude(i,j)=-168.995+(0.01*(i-1))
              lon(i)=-168.995+(0.01*(i-1))
          end do
      end do

! ============================================================
! BUILD 4 KM GRID
! Same as original
! ============================================================
      do i=1,NX1
          do j=1,NY1
              latitude1(i,j)=7.0+(dx1*(j-1))
              longitude1(i,j)=-169+(dx1*(i-1))
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
      call check(nf90_inquire_variable(ncid,totprecip_id,dimids=dimIDs1))
      call check(nf90_inquire_dimension(ncid,dimIDs1(1),len=numLons))
      call check(nf90_inquire_dimension(ncid,dimIDs1(2),len=numLats))
      dimids1=(/numLons,numLats/)
      call check(nf90_get_var(ncid,totprecip_id,mask))
      call check(nf90_close(ncid))

! ============================================================
! YEAR LOOP
! Same structure as original, but year range is user-defined
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
! Same structure as original, but month range is user-defined
! ============================================================
          DO imonth=start_month,end_month

              WRITE(mm,'(I2.2)') imonth
              yyyymm = yyyy//mm

! ============================================================
! READ MONTHLY MODIS CFM 1 KM
! ============================================================
              filename1='CFM_'//yyyymm//'.nc'
              filename2=trim(cfm_modis_dir)//'/'//trim(filename1)

              write(*,*) "Reading MODIS CFM: ", trim(filename2)

              call check(NF90_OPEN(trim(filename2),NF90_NOWRITE,ncid))
              call check(nf90_inquire(ncid,ndimsp,nvarsp,nattsp,unlimdimidp))
              dimid=1
              call check(nf90_inq_varid(ncid,"CFM",totprecip_id))
              call check(nf90_inquire_variable(ncid,totprecip_id,dimids=dimIDs1))
              call check(nf90_inquire_dimension(ncid,dimIDs1(1),len=numLons))
              call check(nf90_inquire_dimension(ncid,dimIDs1(2),len=numLats))
              dimids1=(/numLons,numLats/)
              call check(nf90_get_var(ncid,totprecip_id,cfm_1km_mod))
              call check(nf90_close(ncid))

! ============================================================
! READ MONTHLY MYD CFM 1 KM
! ============================================================
              filename2=trim(cfm_myd_dir)//'/'//trim(filename1)

              write(*,*) "Reading MYD CFM: ", trim(filename2)

              call check(NF90_OPEN(trim(filename2),NF90_NOWRITE,ncid))
              call check(nf90_inquire(ncid,ndimsp,nvarsp,nattsp,unlimdimidp))
              dimid=1
              call check(nf90_inq_varid(ncid,"CFM",totprecip_id))
              call check(nf90_inquire_variable(ncid,totprecip_id,dimids=dimIDs1))
              call check(nf90_inquire_dimension(ncid,dimIDs1(1),len=numLons))
              call check(nf90_inquire_dimension(ncid,dimIDs1(2),len=numLats))
              dimids1=(/numLons,numLats/)
              call check(nf90_get_var(ncid,totprecip_id,cfm_1km_myd))
              call check(nf90_close(ncid))

! ============================================================
! DAY LOOP
! Same as original
! ============================================================
              DO iday=1,nday(imonth)

                  WRITE(dd,'(I2.2)') iday

                  write(*,*) iyear,imonth,iday

! ============================================================
! DAILY FILE NAME
! Same original file naming convention:
! LIS_HIST_YYYYMMDD0000.d01.nc
! ============================================================
                  filename1='LIS_HIST_'//yyyymm//dd//'0000.d01.nc'

! ============================================================
! READ DAILY NLDAS3 4 KM / 10 KM BRASETH FILE
! ============================================================
                  fileplacea=trim(nldas3_input_root)//'/'//yyyymm
                  filename1a=trim(fileplacea)//'/'//trim(filename1)

                  write(*,*) "Reading NLDAS3 input: ", trim(filename1a)

                  call check(NF90_OPEN(trim(filename1a), NF90_NOWRITE, ncid))
                  call check(nf90_inquire(ncid, ndimsp, nvarsp, nattsp, unlimdimidp))
                  dimid=1
                  call check(nf90_inq_varid(ncid, "TotalPrecip_tavg", totprecip_id))
                  call check(nf90_inquire_variable(ncid, totprecip_id, dimids = dimIDs))
                  call check(nf90_inquire_dimension(ncid, dimIDs(1), len = numLons))
                  call check(nf90_inquire_dimension(ncid, dimIDs(2), len = numLats))
                  dimids = (/numLons, numLats/)
                  call check(nf90_get_var(ncid, totprecip_id,nldas3_4km))
                  call check(nf90_close(ncid))

! ============================================================
! READ DAILY MERRA2 4 KM FILE
! ============================================================
                  fileplacea=trim(merra2_input_root)//'/'//yyyymm
                  filename1a=trim(fileplacea)//'/'//trim(filename1)

                  write(*,*) "Reading MERRA2 4km input: ", trim(filename1a)

                  call check(NF90_OPEN(trim(filename1a), NF90_NOWRITE, ncid))
                  call check(nf90_inquire(ncid, ndimsp, nvarsp, nattsp, unlimdimidp))
                  dimid=1
                  call check(nf90_inq_varid(ncid, "TotalPrecip", totprecip_id))
                  call check(nf90_inquire_variable(ncid, totprecip_id, dimids = dimIDs))
                  call check(nf90_inquire_dimension(ncid, dimIDs(1), len = numLons))
                  call check(nf90_inquire_dimension(ncid, dimIDs(2), len = numLats))
                  dimids = (/numLons, numLats/)
                  call check(nf90_get_var(ncid, totprecip_id,merra2_4km))
                  call check(nf90_close(ncid))

! ============================================================
! STEP 1:
! Build 1 km CFM and interpolate/fill NLDAS3 to 1 km
! ============================================================
                  nldas3=0.0d0

                  do i=1,NX
                      do j=1,NY

                          cfm_1km(i,j)=(cfm_1km_myd(i,j)+cfm_1km_mod(i,j))/2.0d0

                          inv_dist_pt=0.0d0

                          xi=int((longitude(i,j)-(-169.0d0))/dx1)+1
                          yi=int((latitude(i,j)-7.0d0)/dx1)+1

                          if (mask(i,j)==1.0d0 .and. nldas3_4km(xi,yi)<0.0d0) then
                              nldas3_4km(xi,yi)=merra2_4km(xi,yi)*86400.0d0
                          end if

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

                                          if (nldas3_4km(x,y)>=0.0d0) then

                                              dist=great_circle_distance(latb,lonb,lata,lona)

                                              if (dist==0.0d0) then
                                                  inv_dist=1.0d0
                                              else
                                                  inv_dist=(1.0d0/dist)**inv_p
                                              end if

                                              inv_dist_pt=inv_dist_pt+inv_dist
                                              nldas3(i,j)=nldas3(i,j)+inv_dist*nldas3_4km(x,y)

                                          end if

                                      end if
                                  end if

                              end do
                          end do

                          if (inv_dist_pt==0.0d0) then
                              nldas3(i,j)=-9999.0d0
                          else
                              nldas3(i,j)=nldas3(i,j)/inv_dist_pt
                          end if

                      end do
                  end do

! ============================================================
! STEP 2:
! Aggregate 1 km NLDAS3 × CFM to 4 km
! ============================================================
                  pcfm_4km=0.0d0

                  do i=1,NX1
                      do j=1,NY1

                          inv_dist_4km=0.0d0

                          xi=int((longitude1(i,j)-(-168.995d0))/dx)+1
                          yi=int((latitude1(i,j)-7.005d0)/dx)+1

                          lata=latitude1(i,j)
                          lona=longitude1(i,j)

                          ii1=xi-5
                          ii2=xi+5
                          jj1=yi-5
                          jj2=yi+5

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

                                  if (cfm_1km(x,y)>=0.0d0 .and. nldas3(x,y)>=0.0d0) then

                                      if (abs(latitude1(i,j)-latitude(x,y))<=dx1) then
                                          if (abs(abs(longitude1(i,j))-abs(longitude(x,y)))<=dx1) then

                                              if (nldas3(x,y)>=0.0d0 .and. cfm_1km(x,y)>0.0d0) then

                                                  latb=latitude(x,y)
                                                  lonb=longitude(x,y)

                                                  dist=great_circle_distance(latb,lonb,lata,lona)

                                                  if (dist==0.0d0) then
                                                      inv_dist=1.0d0
                                                  else
                                                      inv_dist=(1.0d0/dist)**inv_p
                                                  end if

                                                  inv_dist_4km=inv_dist_4km+inv_dist
                                                  pcfm_4km(i,j)=pcfm_4km(i,j)+inv_dist*nldas3(x,y)*cfm_1km(x,y)

                                              end if

                                          end if
                                      end if

                                  end if

                              end do
                          end do

                          if (inv_dist_4km==0.0d0) then
                              pcfm_4km(i,j)=0.0d0
                          else
                              pcfm_4km(i,j)=pcfm_4km(i,j)/inv_dist_4km
                          end if

                      end do
                  end do

! ============================================================
! STEP 3:
! Apply cloud-frequency downscaling at 1 km
! ============================================================
                  do i=1,NX
                      do j=1,NY

                          if (cfm_1km(i,j)>=0.0d0 .and. nldas3(i,j)>=0.0d0) then

                              inv_dist_pt=0.0d0
                              nldas3_pt=0.0d0
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

                                              if (nldas3_4km(x,y)>=0.0d0) then

                                                  dist=great_circle_distance(latb,lonb,lata,lona)

                                                  if (dist==0.0d0) then
                                                      inv_dist=1.0d0
                                                  else
                                                      inv_dist=(1.0d0/dist)**inv_p
                                                  end if

                                                  inv_dist_pt=inv_dist_pt+inv_dist
                                                  nldas3_pt=nldas3_pt+inv_dist*nldas3_4km(x,y)
                                                  cfm_pt=cfm_pt+inv_dist*pcfm_4km(x,y)

                                              end if

                                          end if
                                      end if

                                  end do
                              end do

                              if (inv_dist_pt==0.0d0) then
                                  nldas3_pt=0.0d0
                                  cfm_pt=0.0d0
                              else
                                  nldas3_pt=nldas3_pt/inv_dist_pt
                                  cfm_pt=cfm_pt/inv_dist_pt
                              end if

                              if (cfm_pt==0.0d0) then
                                  nldas3_hr=nldas3(i,j)
                              else
                                  nldas3_hr=((cfm_1km(i,j)*nldas3(i,j))/cfm_pt)*nldas3_pt
                              end if

                              if (nldas3_hr>(1.5d0*nldas3_pt)) nldas3_hr=nldas3(i,j)

                          else

                              nldas3_hr=nldas3(i,j)

                          end if

                          nldas3_1km(i,j)=nldas3_hr

                      end do
                  end do

! ============================================================
! STEP 4:
! ============================================================
                  do i=1,NX
                      do j=1,NY
                          if (mask(i,j)==1.0d0 .and. nldas3_1km(i,j)<0.0d0) then
                              stop
                          end if
                      end do
                  end do

! ============================================================
! WRITE DAILY 1 KM OUTPUT
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

                  dimids1=(/lon_dimid,lat_dimid/)

                  call check(nf90_def_var(ncid,"lat",NF90_REAL,dimids1,latitude_id))
                  CALL check(nf90_put_att(ncid,latitude_id,"units","degree_north"))
                  call check(nf90_put_att(ncid,latitude_id,"_FillValue", -9999.00))
                  CALL check(nf90_put_att(ncid,latitude_id,"long_name","latitude"))

                  call check(nf90_def_var(ncid,"lon",NF90_REAL,dimids1,longitude_id))
                  CALL check(nf90_put_att(ncid,longitude_id,"units","degree_east"))
                  call check(nf90_put_att(ncid,longitude_id,"_FillValue", -9999.00))
                  CALL check(nf90_put_att(ncid,longitude_id,"long_name","longitude"))

                  dimids=(/lon_dimid,lat_dimid/)

                  call check(nf90_def_var(ncid,"TotalPrecip_tavg",NF90_REAL,dimids,totprecip_id))
                  CALL check(nf90_put_att(ncid,totprecip_id,"units","mm"))
                  call check(nf90_put_att(ncid,totprecip_id,"_FillValue", -9999.00))
                  CALL check(nf90_put_att(ncid,totprecip_id,"long_name","total precipitation"))

                  CALL check(nf90_enddef(ncid))

                  call check(nf90_put_var(ncid,lat_id,lat))
                  call check(nf90_put_var(ncid,lon_id,lon))
                  call check(nf90_put_var(ncid,totprecip_id,nldas3_1km))
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

end Program downscaled_nldas3_cloud_1km
