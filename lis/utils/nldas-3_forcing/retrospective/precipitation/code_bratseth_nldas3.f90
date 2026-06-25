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
! !ROUTINE: Pcp_assimilation
!
! !REVISION HISTORY:
!  12 Jun 2026  Fadji Maina; Initial specification
!
! !COMPILATION:
!  This program can be compiled on NASA Discover using the Intel MPI Fortran
!  compiler and the LISF NetCDF/HDF5 libraries as follows:
!
!  mpif90 -g -check all -traceback -names lowercase -convert big_endian 
!    -assume byterecl 
!    -I/discover/nobackup/projects/lis/libs/sles-12.3/netcdf/4.8.1_intel-2021.4.0/include 
!    code_bratseth_nldas3.f90 -o code_bratseth_nldas3
!    -L/discover/nobackup/projects/lis/libs/sles-12.3/netcdf/4.8.1_intel-2021.4.0/lib 
!    -L/discover/nobackup/projects/lis/libs/sles-12.3/hdf5/1.12.1_intel-2021.4.0/lib 
!    -lnetcdff -lnetcdf -lhdf5_hl -lhdf5 -lz -lcurl 
!    -Wl,--no-relax -shared-intel
!
! !DESCRIPTION:
!  This program assimilates IMERG and CAPA precipitation estimates into a
!  MERRA-2/LIS precipitation background using an MPI-parallel BRASH analysis.
!  The BRASH algorithm is based on the Bratseth analysis approach and follows
!  the precipitation merging framework described by Kemp et al. (2022). The
!  output is a daily NetCDF precipitation analysis on the NLDAS-3 grid.
!
!EOP

Program Pcp_assimilation
use mpi
use netcdf
  implicit none
      !include 'mpif.h'
      integer i,j,n,iyear,imonth,iday,nday(12)
      integer comm, myid, numprocs,ierr,npts,npts_mpi,imaxabsdiff,ip,nobs,jd,id,xi,yi
      integer, parameter:: NX=2926, NY=1626,NX2=187,NY2=111,NX3=118,NY3=66
      Character (LEN=200) filename,fileplace,fvhour2,fvhour1,date_str,iav
      Character (LEN=200) filename1,fileplace1,fileplace2,fvario,filenamenc
      Integer, parameter :: ndims=2,ndims1=2,NX4=11700,NY4=6500,NX1=1171,NY1=651
      integer :: dimids(ndims),numLons,numLats,dimid,dimids1(ndims1)
      real*8 ::  merra2ErrScaleLength,convergeThresh,sigmasqr_back,sigmasqr_obs1
      real*8 :: sigmasqr_obs2,backErrScaleLength,maxabsdiff,y_prev,y_new
      real*8 :: sigmasqr_merra2,sigmasqr_obscapa,sigmasqr_obsimerg,sigmasqr_imerg,sigmasqr_capa,sigmasqr_obs
      real*8 :: obs1ErrScaleLength,obs2ErrScaleLength,sigmasqr_satimerg,sigmasqr_satcapa,imergErrScaleLength,capaErrScaleLength
      real*8, allocatable, dimension (:,:) :: latimerg,lonimerg,latitude,longitude,precepimerg,precepmerra2,precep,lat_1km,lon_1km
      real*8, allocatable, dimension (:,:) :: obs_typxy,apcpobs,apcpback,mask,buffer,ele,mask_1km,precepcapa,precepcfia
      real*8, allocatable, dimension (:) :: lat1,lon1,precep_back,precep_obs,precep_obs1,precep_obs2,pprev_est,pprev_ana,pnew_est_mpi,pnew_ana_mpi
      real*8, allocatable, dimension (:) :: pnew_est,pnew_ana,sumObsEstimates_mpi
      real*8, allocatable, dimension (:) :: invDataDensities,sumObsEstimates,mrgp_pcp,mrgp_pcp_mpi,invDataDensities_mpi
      real*8, allocatable, dimension (:) :: lat,lon,latobs,lonobs
      integer, allocatable, dimension (:,:) :: id_imerg,jd_imerg
      integer :: ncid,totprecip_id,lat_dimid,lon_dimid,lat_id,lon_id,latitude_id,longitude_id,length_id,sigma_id,obs_id,apcpobs_id,apcpback_id
      integer :: x, y, ndimsp, nvarsp, nattsp, unlimdimidp,npasses,corrobs,ccount,ele_id
      integer, allocatable, dimension (:) :: start,count,id_capa,jd_capa,obs_typ
      logical :: done
	npts=nx*ny

	allocate (start(ndims),count(ndims),precepimerg(NX1,NY1),precepmerra2(NX,NY),precep(NX,NY),obs_typ(npts),obs_typxy(NX,NY))
        allocate (precepcapa(NX1,NY1),precepcfia(NX1,NY1))
        allocate (latitude(NX,NY),longitude(NX,NY),lat(ny),lon(nx),apcpobs(NX,NY))
	allocate (apcpback(NX,NY),lat1(npts),lon1(npts),id_imerg(NX1,NY1),jd_imerg(NX1,NY1),latimerg(NX1,NY1),lonimerg(NX1,NY1))
	allocate (precep_back(npts),pprev_est(npts),pprev_ana(npts),sumObsEstimates_mpi(npts))
        allocate (invDataDensities(npts),sumObsEstimates(npts),mrgp_pcp(npts),invDataDensities_mpi(npts),ele(nx3,ny3))
        allocate (lat_1km(NX4,NY4),lon_1km(NX4,NY4),mask(NX,NY),precep_obs1(npts),precep_obs2(npts))
        allocate (mrgp_pcp_mpi(npts),pnew_est_mpi(npts),pnew_ana_mpi(npts),pnew_est(npts),pnew_ana(npts),buffer(nx3,ny3),mask_1km(NX4,NY4))

      fileplace='2000'

      nday = (/31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/)

      nobs=npts
      allocate (precep_obs(nobs),latobs(nobs),lonobs(nobs))

        call MPI_INIT(ierr)
        call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
        call MPI_COMM_SIZE( MPI_COMM_WORLD,numprocs,ierr)


! Read LANDMASK from 4km domain
filename = 'nldas_domain_4km_0823.nc'
call check(NF90_OPEN(filename, NF90_NOWRITE, ncid))
call check(nf90_inquire(ncid, ndimsp, nvarsp, nattsp, unlimdimidp))
call check(nf90_inq_varid(ncid, "LANDMASK", totprecip_id))
call check(nf90_inquire_variable(ncid, totprecip_id, dimids=dimIDs))
call check(nf90_inquire_dimension(ncid, dimIDs(1), len=numLons))
call check(nf90_inquire_dimension(ncid, dimIDs(2), len=numLats))
call check(nf90_get_var(ncid, totprecip_id, mask))
call check(nf90_close(ncid))

! Read LANDMASK from 1km domain
filename = 'lis_input.nldas3.noahmp401.1km.nc'
call check(NF90_OPEN(filename, NF90_NOWRITE, ncid))
call check(nf90_inquire(ncid, ndimsp, nvarsp, nattsp, unlimdimidp))
call check(nf90_inq_varid(ncid, "LANDMASK", totprecip_id))
call check(nf90_inquire_variable(ncid, totprecip_id, dimids=dimIDs1))
call check(nf90_inquire_dimension(ncid, dimIDs1(1), len=numLons))
call check(nf90_inquire_dimension(ncid, dimIDs1(2), len=numLats))
call check(nf90_get_var(ncid, totprecip_id, mask_1km))
call check(nf90_close(ncid))

! Read LANDMASK and ELEVATION from 100km domain
filename = 'nldas_domain_100km_0823.nc'
call check(NF90_OPEN(filename, NF90_NOWRITE, ncid))
call check(nf90_inquire(ncid, ndimsp, nvarsp, nattsp, unlimdimidp))
call check(nf90_inq_varid(ncid, "LANDMASK", totprecip_id))
call check(nf90_inquire_variable(ncid, totprecip_id, dimids=dimIDs))
call check(nf90_inquire_dimension(ncid, dimIDs(1), len=numLons))
call check(nf90_inquire_dimension(ncid, dimIDs(2), len=numLats))
call check(nf90_inq_varid(ncid, "ELEVATION", ele_id))
call check(nf90_get_var(ncid, totprecip_id, buffer))
call check(nf90_get_var(ncid, ele_id, ele))
call check(nf90_close(ncid))

        do i=1,nx4
                do j=1,ny4
                        lat_1km(i,j)=7.005+(0.01*(j-1))
                        lon_1km(i,j)=-168.995+(0.01*(i-1))
                        xi=int((lon_1km(i,j)-(-169))/0.04)+1
                        yi=int((lat_1km(i,j)-7.0)/0.04)+1
                        if (mask_1km(i,j)==1 .and. mask(xi,yi)==0) then
                                mask(xi,yi)=1
                        end if
                end do
        end do


        do i=1,nx
                do j=1,ny
                        latitude(i,j)=7.0+(0.04*(j-1))
                        lat(j)=7.0+(0.04*(j-1))
                        longitude(i,j)=-169+(0.04*(i-1))
                        lon(i)=-169+(0.04*(i-1))
                        if (int(mask(i,j))==0) then
                                if (i>1073 .and. i<1311) then
                                        if (j>1428 .and. j<1518) then
                                                mask(i,j)=1.0
                                        end if
                                end if
                                if (i>1275 .and. i<1833) then
                                        if (j>1020 .and. j<1430) then
                                                mask(i,j)=1.0
                                        end if
                                end if
                                if (i>1903 .and. i<2345) then
                                        if (j>830 .and. j<1094) then
                                                mask(i,j)=1.0
                                        end if
                                end if
                                if (i>2578 .and. i<2653) then
                                        if (j>1147 .and. j<1195) then
                                                mask(i,j)=1.0
                                        end if
                                end if
                                id=int((longitude(i,j)-(-169))/1.0)+1
                                jd=int((latitude(i,j)-7.0)/1.0)+1

                                if (int(mask(i,j))==0 .and. ele(id,jd) .ne. -9999.0) then
                                        mask(i,j)=1.0
                                end if
                        end if
                end do
        end do

do iyear = 26, 26
    ! Adjust February days for leap year
    if (mod(2000 + iyear, 4) == 0) then
        nday(2) = 29
    else
        nday(2) = 28
    end if

    do imonth = 3, 3
        do iday = 27, nday(imonth)
            call MPI_BARRIER(MPI_COMM_WORLD, ierr)

            ! Compose date string YYYYMMDDHH with HH='00'
            write(date_str, '(I4.4,I2.2,I2.2,A2)') 2000 + iyear, imonth, iday, '00'

            filename = 'LIS_HIST_' // trim(date_str) // '00.d01.nc'       !trim(date_str(1:8))// '00.d01.nc'
            fileplace=trim(date_str(1:6))
!            if (myid == 0) write(*,*) filename 
            ! Read statistics files
            fvario = "./merra2_semivar/fit_variomerr_" // trim(date_str(1:8)) // ".txt"
            open(11, file=trim(fvario), status='old', action='read', iostat=ierr)
                if (ierr /= 0) then
                 if (myid == 0) write(*,*) 'Error opening file: ', trim(fvario)
                 stop
                endif
            read(11, *) iav, sigmasqr_obs
            read(11, *) iav, sigmasqr_merra2
            read(11, *) iav, merra2ErrScaleLength
            close(11)

            fvario = "./gimerg10km_semivar/fit_varioimer_" // trim(date_str(1:8)) // ".txt"
            open(11, file=trim(fvario), status='old', action='read', iostat=ierr)
                if (ierr /= 0) then
                 if (myid == 0) write(*,*) 'Error opening file: ', trim(fvario)
                 stop
                endif
            read(11, *) iav, sigmasqr_obsimerg
            read(11, *) iav, sigmasqr_satimerg
            read(11, *) iav, ImergErrScaleLength
            close(11)

            fvario = "./capa10km_semivar/fit_variocapa_" // trim(date_str(1:8)) // ".txt"
            open(11, file=trim(fvario), status='old', action='read', iostat=ierr)
                if (ierr /= 0) then
                 if (myid == 0) write(*,*) 'Error opening file: ', trim(fvario)
                 stop
                endif
            read(11, *) iav, sigmasqr_obscapa
            read(11, *) iav, sigmasqr_satcapa
            read(11, *) iav, capaErrScaleLength
            close(11)

            ! Compute adjusted variances
            sigmasqr_capa = sigmasqr_obs * (sigmasqr_satcapa / sigmasqr_obscapa)
            sigmasqr_imerg = sigmasqr_obs * (sigmasqr_satimerg / sigmasqr_obsimerg)

            ! Define helper subroutine to open NetCDF, get dimensions and variable

            call read_netcdf("./merra2_daily_4km/" // trim(adjustl(fileplace)) // "/",  trim(adjustl(filename)), &
                             "TotalPrecip", precepmerra2, ncid, dimIDs, numLons, numLats)

            call read_netcdf("./gimerg10km_daily/" // trim(adjustl(fileplace)) // "/", filename, &
                             "TotalPrecip", precepimerg, ncid, dimIDs, numLons, numLats)

            call read_netcdf("./capar10km_daily/" // trim(adjustl(fileplace)) // "/", filename, &
                             "TotalPrecip_tavg", precepcapa, ncid, dimIDs, numLons, numLats)

            call read_netcdf("./cfia10km_daily/" // trim(adjustl(fileplace)) // "/", filename, &
                             "TotalPrecip_tavg", precepcfia, ncid, dimIDs, numLons, numLats)

            ! Initialize precipitation arrays and coordinates
            precep_back = -9999.0
            precep_obs = -9999.0

                        do i=1,nx1
                                do j=1,ny1
                                        latimerg(i,j)=7+(j-1)*0.1
                                        lonimerg(i,j)=-169+(i-1)*0.1
                                        id_imerg(i,j)=int((lonimerg(i,j)-(-169))/0.04)+1
                                        jd_imerg(i,j)=int((latimerg(i,j)-7.0)/0.04)+1
                                end do
                        end do

            if (myid == 0) write(*,*) iyear, imonth, iday
            ip = 0
            do i = 1, nx
                do j = 1, ny
                    ip = ip + 1
                    if (int(mask(i,j)) /= 0) then
                        precep_back(ip) = precepmerra2(i,j) * 86400.0
                    else
                        precep_back(ip) = -9999.0
                    end if
                    pprev_ana(ip) = precep_back(ip)
                    pprev_est(ip) = precep_back(ip)
                    lat1(ip) = latitude(i,j)
                    lon1(ip) = longitude(i,j)
                end do
            end do

                        obs_typ=0
                        precep_obs1=-9999.0  ;   precep_obs2=-9999.0  ; latobs=0.0  ; lonobs=0.0
                        do i=1,nx1
                                do j=1,ny1
                                        ip=((id_imerg(i,j)-1)*ny)+jd_imerg(i,j)
                                        latobs(ip)=latimerg(i,j)
                                        lonobs(ip)=lonimerg(i,j)
                                        precep_obs1(ip)=-9999.0
                                        precep_obs2(ip) = -9999.0
                                        if (precepimerg(i,j)<0.0) then
                                                precep_obs1(ip)=-9999.0
                                        else
                                                precep_obs1(ip)=precepimerg(i,j)*86400
                                        end if
                                        if (precepcfia(i,j) > 0.7) then
                                                precep_obs2(ip) = precepcapa(i,j)*86400.0
                                                if (precep_obs2(ip)<0.0) precep_obs2(ip)=-9999.0
                                        else
                                                precep_obs2(ip) = -9999.0
                                        end if
                                end do
                        end do

            ! Set variance and error scale length variables
            sigmasqr_back = sigmasqr_merra2
            backErrScaleLength = merra2ErrScaleLength * 1000.0
            sigmasqr_obs1 = sigmasqr_imerg
            sigmasqr_obs2 = sigmasqr_capa
            obs1ErrScaleLength = ImergErrScaleLength * 1000.0
            obs2ErrScaleLength = capaErrScaleLength * 1000.0

            ! Broadcast parameters
            call MPI_BCAST(obs1ErrScaleLength, 1, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(obs2ErrScaleLength, 1, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(backErrScaleLength, 1, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(sigmasqr_back, 1, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(sigmasqr_obs1, 1, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(sigmasqr_obs2, 1, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(npts, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(lat1, npts, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(lon1, npts, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(latobs, npts, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(lonobs, npts, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(precep_obs1, npts, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(precep_obs2, npts, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(obs_typ, npts, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(precep_back, npts, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)

            call MPI_BARRIER(MPI_COMM_WORLD, ierr)

            invDataDensities = 0.0
            invDataDensities_mpi = 0.0

            if (myid == 0) write(*,*) 'call calc_invDataDensities'
!            if (myid == 0) write(*,*) sigmasqr_back, sigmasqr_obs1, sigmasqr_obs2,backErrScaleLength, obs1ErrScaleLength, obs2ErrScaleLength
            call calc_invDataDensities(myid, numprocs, npts, obs_typ, sigmasqr_back, sigmasqr_obs1, sigmasqr_obs2, &
                                       backErrScaleLength, obs1ErrScaleLength, obs2ErrScaleLength, lat1, lon1, latobs, lonobs, &
                                       invDataDensities_mpi, precep_obs1,precep_obs2, precep_back)


            invDataDensities = 0.0
            call MPI_BARRIER(MPI_COMM_WORLD, ierr)
            call MPI_ALLREDUCE(invDataDensities_mpi, invDataDensities, npts, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)

            ! Iterative observation analysis
            convergeThresh = 0.1
            if (myid == 0) write(*,*) 'call calc_obsAnalysis'
            npasses = 0
            sumObsEstimates = 0.0
            sumObsEstimates_mpi = 0.0

            do
                call MPI_BARRIER(MPI_COMM_WORLD, ierr)
                call calc_obsAnalysis(myid, numprocs, npts, obs_typ, sigmasqr_obs1, sigmasqr_obs2, sigmasqr_back, invDataDensities, &
                                      backErrScaleLength, obs1ErrScaleLength, obs2ErrScaleLength, npasses, pprev_ana, pprev_est, &
                                      pnew_ana_mpi, pnew_est_mpi, sumObsEstimates_mpi, lat1, lon1, latobs, lonobs, precep_obs1,precep_obs2, precep_back)

                pnew_est = 0.0
                pnew_ana = 0.0
                sumObsEstimates = 0.0

                call MPI_BARRIER(MPI_COMM_WORLD, ierr)
                call MPI_ALLREDUCE(pnew_ana_mpi, pnew_ana, npts, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
                call MPI_ALLREDUCE(pnew_est_mpi, pnew_est, npts, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
                call MPI_ALLREDUCE(sumObsEstimates_mpi, sumObsEstimates, npts, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)

                done = .true.
                maxabsdiff = 0.0
                imaxabsdiff = 0
                npasses = npasses + 1

                call MPI_BARRIER(MPI_COMM_WORLD, ierr)

                do i = 1, npts
                    if (precep_back(i) /= -9999.0) then
                        y_prev = pprev_ana(i)
                        y_new = pnew_ana(i)
                        if (abs(y_prev - y_new) > convergeThresh) then
                            if (abs(y_prev - y_new) > maxabsdiff) then
                                maxabsdiff = abs(y_prev - y_new)
                                imaxabsdiff = i
                            end if
                            done = .false.
                        end if
                        pprev_est(i) = pnew_est(i)
                        pprev_ana(i) = pnew_ana(i)
                    end if
                end do

                if (myid == 0) then
                    write(*,*) 'npasses', npasses, maxabsdiff, imaxabsdiff
                end if

                if (npasses > 5 .or. done) exit

                call MPI_BARRIER(MPI_COMM_WORLD, ierr)
            end do

            call MPI_BARRIER(MPI_COMM_WORLD, ierr)
            call MPI_BCAST(npasses, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            call MPI_BARRIER(MPI_COMM_WORLD, ierr)

            if (myid == 0) write(*,*) 'call calc_gridAnalysis'
            mrgp_pcp_mpi = 0.0
            call calc_gridAnalysis(myid, numprocs, npts, obs_typ, sigmasqr_back, sigmasqr_obs1, sigmasqr_obs2, backErrScaleLength, obs1ErrScaleLength, &
                                   obs2ErrScaleLength, precep_obs1, precep_obs2, precep_back, sumObsEstimates, npasses, invDataDensities, lat1, lon1, latobs, lonobs, mrgp_pcp_mpi)

            mrgp_pcp = 0.0
            call MPI_BARRIER(MPI_COMM_WORLD, ierr)
            call MPI_ALLREDUCE(mrgp_pcp_mpi, mrgp_pcp, npts, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)

            if (myid == 0) then
                ip = 0
                apcpobs = 0.0
                precep = 0.0
                obs_typxy = 0.0
                apcpback = 0.0

                do i = 1, nx
                    do j = 1, ny
                        ip = ip + 1
                        precep(i,j) = mrgp_pcp(ip)
                        if (mrgp_pcp(ip)>1.0 .and. precep_back(ip)>1.0) then
                                apcpback(i,j) = mrgp_pcp(ip)/precep_back(ip)
                        else
                                apcpback(i,j) =0.0
                        end if
                        if (precep(i,j) < 0.0) then
                                precep(i,j) = -9999.0
                                apcpback(i,j) =-9999.0
                        end if
                        if (precep_back(ip) < 0.0) then
                                precep(i,j) = -9999.0
                                apcpback(i,j) =-9999.0
                        end if
                    end do
                end do

                ! Create output directory if it does not exist
                fileplace1 = "./nldas3_10km_braseth/" // trim(fileplace) // "/"
                call execute_command_line("mkdir -p " // trim(fileplace1))

                filenamenc = trim(adjustl(fileplace1)) // "/" // filename

                ! Create NetCDF file and define variables (use previously defined subroutine for clarity)
                call create_and_write_netcdf(filenamenc, lat, lon, precep, apcpback, latitude, longitude)
            end if

            call MPI_BARRIER(MPI_COMM_WORLD, ierr)

        end do ! iday
    end do ! imonth
end do ! iyear


14 FORMAT(2000(1x,e13.7))


contains
  subroutine check(status)
    integer, intent ( in) :: status

    if(status /= nf90_noerr) then
      print *, trim(nf90_strerror(status))
      stop "Stopped"
    end if
  end subroutine check


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
   subroutine calc_invDataDensities(myid,numprocs,npts,obs_typ,sigmasqr_back,sigmasqr_obs1,sigmasqr_obs2, &
        backErrScaleLength,obs1ErrScaleLength,obs2ErrScaleLength,lat,lon,latobs,lonobs,invDataDensities,precep_obs1,precep_obs2,precep_back)
      ! Arguments
      real*8, allocatable, intent(in) :: lat(:),lon(:),latobs(:),lonobs(:)
      integer, allocatable, intent(in) :: obs_typ(:)
      integer, intent(in) :: npts,myid, numprocs
      real*8, intent(in) :: sigmasqr_back,sigmasqr_obs1,sigmasqr_obs2,backErrScaleLength,obs1ErrScaleLength,obs2ErrScaleLength
      real*8, allocatable, intent(in) :: precep_back(:),precep_obs1(:),precep_obs2(:)
      real*8, allocatable :: invDataDensities(:)
      ! Local variables
      real*8 :: dist,lat1,lon1,lat2_1,lon2_1,lat2_2,lon2_2,sigmaOSqriob,dust,o1,max_dist,data_sum,sigmaOSqri_1,sigmaOSqri_2
      real*8 :: b, num, denom,sigmaOSqrj1,sigmaOSqri,sigmaBSqr,bErrScaleLength,inv_scale_length,oErrScaleLength
      integer :: i,j,id,id_incr,i1,j1,j2,xi,yi,xj,yj,jj1,jj2,ii1,ii2,x,y,ic
      real*8,allocatable, dimension (:)  ::dataDensities
      integer, parameter:: NX=2926, NY=1626

allocate(dataDensities(npts))
dataDensities = 0.0
id = -1
id_incr = 1

sigmaBSqr = sigmasqr_back
bErrScaleLength = backErrScaleLength

if (obs1ErrScaleLength == 0.0) then
  max_dist = 2 * MAX(backErrScaleLength, obs2ErrScaleLength) / 1000.0
else if (obs2ErrScaleLength == 0.0) then
  max_dist = 2 * MAX(backErrScaleLength, obs1ErrScaleLength) / 1000.0
else
  max_dist = 2 * MAX(backErrScaleLength, obs1ErrScaleLength, obs2ErrScaleLength) / 1000.0
end if

ic = int(max_dist / 4.0)

do i = 1, npts
if (precep_back(i) < 0.0 .or. precep_obs1(i)<0.0) cycle  ! skip invalid points
  call pick_proc(id, id_incr, numprocs)
  if (id /= myid) cycle
  if (precep_back(i) >= 0.0 .and. precep_obs1(i)>=0.0) then
    ! For obs1 (e.g. IMERG) at this grid point
    if (precep_obs1(i) >= 0.0) then
      sigmaOSqri_1 = sigmasqr_obs1
      lat2_1 = latobs(i)
      lon2_1 = lonobs(i)
    else
      sigmaOSqri_1 = -1.0   ! Flag invalid, won't be used
    end if

    ! For obs2 (e.g. CAPA) at this grid point with small lon offset
    if (precep_obs2(i) >= 0.0) then
      sigmaOSqri_2 = sigmasqr_obs2
      lat2_2 = latobs(i) - 0.02d0
      lon2_2 = lonobs(i) - 0.02d0  ! small longitude offset for obs2
    else
      sigmaOSqri_2 = -1.0
    end if

    xi = int((lon(i) - (-169.0d0)) / 0.04d0) + 1
    yi = int((lat(i) - 7.0d0) / 0.04d0) + 1
    ii1 = max(1, xi - ic)
    ii2 = min(nx, xi + ic)
    jj1 = max(1, yi - ic)
    jj2 = min(ny, yi + ic)

    do x = ii1, ii2
      do y = jj1, jj2
        j1 = ((x - 1) * ny) + y

        if (precep_back(j1) >= 0.0) then
          ! Contribution from obs1 at j1 (no offset)
          if (precep_obs1(j1) >= 0.0) then
            lat1 = latobs(j1)
            lon1 = lonobs(j1)
            sigmaOSqrj1 = sigmasqr_obs1
            oErrScaleLength = obs1ErrScaleLength
            dist = 0.0d0
            if ((precep_obs1(i) >= 0.0) .and. (precep_obs1(j1) >= 0.0)) then
              if (i /= j1) then
                dist = great_circle_distance(lat2_1, lon2_1, lat1, lon1)
              end if
              if (dist * 0.001d0 < max_dist) then
                b = backErrCov(sigmaBSqr, dist, bErrScaleLength)
                num = b
                if (i == j1) then
                  num = num + sigmaOSqrj1
                else
                  o1 = obsErrCov(sigmaOSqrj1, oErrScaleLength, dist)
                  num = num + o1
                end if
                denom = sqrt((sigmaBSqr + sigmaOSqri_1) * (sigmaBSqr + sigmaOSqrj1))
                dataDensities(i) = dataDensities(i) + (num / denom)
              end if
            end if
          end if

          ! Contribution from obs2 at j1 (with lon offset)
          if (precep_obs2(j1) >= 0.0) then
            lat1 = latobs(j1)- 0.02d0
            lon1 = lonobs(j1) - 0.02d0
            sigmaOSqrj1 = sigmasqr_obs2
            oErrScaleLength = obs2ErrScaleLength

            dist = 0.0d0
            if ((precep_obs2(i) >= 0.0) .and. (precep_obs2(j1) >= 0.0)) then
              if (i /= j1) then
                dist = great_circle_distance(lat2_2, lon2_2, lat1, lon1)
              end if
              if (dist * 0.001d0 < max_dist) then
                b = backErrCov(sigmaBSqr, dist, bErrScaleLength)
                num = b
                if (i == j1) then
                  num = num + sigmaOSqrj1
                else
                  o1 = obsErrCov(sigmaOSqrj1, oErrScaleLength, dist)
                  num = num + o1
                end if
                denom = sqrt((sigmaBSqr + sigmaOSqri_2) * (sigmaBSqr + sigmaOSqrj1))
                dataDensities(i) = dataDensities(i) + (num / denom)
              end if
            end if
          end if

        end if
      end do
    end do

    ! Inverse data density for obs1 and obs2 combined:
    data_sum = 0.0d0
    if (precep_obs1(i) >= 0.0) data_sum = data_sum + (sigmaBSqr + sigmasqr_obs1)
    if (precep_obs2(i) >= 0.0) data_sum = data_sum + (sigmaBSqr + sigmasqr_obs2)

    if (data_sum > 0.0d0) then
      invDataDensities(i) = dataDensities(i) * data_sum
      invDataDensities(i) = 1.0d0 / invDataDensities(i)
    else
      invDataDensities(i) = 0.0d0
    end if

  else
    invDataDensities(i) = 0.0d0
  end if
end do

deallocate(dataDensities)

   end subroutine calc_invDataDensities
   !---------------------------------------------------------------------------
   ! Perform Bratseth analysis at observation points.  Multiple iterations
   ! are made until convergence is reached.  Along the way, the observation
   ! estimates from each iteration are summed for later interpolation to the
   ! grid points.  Note that the *analysis* is also run at the observation
   ! points because in practice the analysis converges before the iterative
   ! "observation estimates" do (see Sashegyi et al 1993).
   ! The observed, background, and analysis values at the observation
   ! points are also collected in an OBA structure for post-processing.
   subroutine calc_obsAnalysis(myid, numprocs,npts,obs_typ,sigmasqr_obs1,sigmasqr_obs2,sigmasqr_back,invDataDensities,&
        backErrScaleLength,obs1ErrScaleLength,obs2ErrScaleLength,npasses,pprev_ana,pprev_est,pnew_ana,pnew_est,&
        sumObsEstimates,lat,lon,latobs,lonobs,precep_obs1,precep_obs2,precep_back)
      implicit none
      ! Arguments
      real*8, allocatable, intent(in) :: lat(:),lon(:),latobs(:),lonobs(:),invDataDensities(:)
      real*8, allocatable, intent(in) :: precep_back(:),precep_obs1(:),precep_obs2(:)
      integer, allocatable, intent(in) :: obs_typ(:)
      real*8, allocatable, intent(inout) :: pnew_ana(:),pprev_ana(:),pnew_est(:),pprev_est(:)
      integer, intent(in) :: npts,npasses,myid, numprocs
      real*8, intent(in) :: sigmasqr_back,sigmasqr_obs1,sigmasqr_obs2,backErrScaleLength,obs1ErrScaleLength,obs2ErrScaleLength
      real*8, allocatable :: sumObsEstimates(:)
      ! Local variables
      real*8 :: dist,lat1,lon1,sigmaOSqri,sigmaOSqrj1,sigmaBSqr,oErrScaleLength,bErrScaleLength,sigmaOSqri_1,sigmaOSqri_2
      real*8 :: b,weight,o1,max_dist,lat2_1,lon2_1,lat2_2,lon2_2
      integer :: i,j,iob,x,y,imaxabsdiff,id,id_incr,i1,j1,xi,yi,xj,yj,jj1,jj2,ii1,ii2,ic
      integer, parameter:: NX=2926, NY=1626
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

	sigmaBSqr=sigmasqr_back
      bErrScaleLength=backErrScaleLength
      id = -1	;	id_incr = 1
      pnew_est=0.0    ;       pnew_ana=0.0
        if (obs1ErrScaleLength==0.0) then
                max_dist=2*(MAX(backErrScaleLength,obs2ErrScaleLength))/1000.0
        else if (obs2ErrScaleLength==0.0) then
                max_dist=2*(MAX(backErrScaleLength,obs1ErrScaleLength))/1000.0
        else
                max_dist=2*(MAX(backErrScaleLength,obs1ErrScaleLength,obs2ErrScaleLength))/1000.0
        end if
      ic=int(max_dist/4.0)
      do i=1,npts
      if (precep_back(i) < 0.0  .or. precep_obs1(i)<0.0) cycle  ! skip invalid points
      call pick_proc(id,id_incr, numprocs)
      if (id .ne. myid) cycle
!      if(precep_back(i)>=0.0) then
       if (precep_back(i) >= 0.0 .and. precep_obs1(i)>=0.0) then
              !For obs1 (e.g. IMERG) at this grid point
        if (precep_obs1(i) >= 0.0) then
                sigmaOSqri_1 = sigmasqr_obs1
                lat2_1 = latobs(i)
                lon2_1 = lonobs(i)
        else
                sigmaOSqri_1 = -1.0   ! Flag invalid, won't be used
        end if

            ! For obs2 (e.g. CAPA) at this grid point with small lon offset
        if (precep_obs2(i) >= 0.0) then
                sigmaOSqri_2 = sigmasqr_obs2
                lat2_2 = latobs(i) - 0.02d0
                lon2_2 = lonobs(i) - 0.02d0  ! small longitude offset for obs2
        else
                sigmaOSqri_1 = -1.0
        end if
                xi=int((lon(i)-(-169))/0.04)+1
                yi=int((lat(i)-7.0)/0.04)+1
                ii1=xi-ic       ;       ii2=xi+ic
                jj1=yi-ic       ;       jj2=yi+ic
                if (ii1<1) ii1=1        ;       if (ii1>nx) ii1=nx
                if (ii2>nx) ii2=nx      ;       if (ii2<1) ii2=1
                if (jj1<1) jj1=1        ;       if (jj1>ny) jj1=ny
                if (jj2>ny) jj2=ny      ;       if (jj2<1) jj2=1
        do x = ii1, ii2
                do y = jj1, jj2
                        j1 = ((x-1)*ny) + y
                        if (precep_back(j1) >= 0.0) then
                        ! --- Obs1 processing ---
                      if (precep_obs1(j1) >= 0.0) then
                        sigmaOSqrj1 = sigmasqr_obs1
                        oErrScaleLength = obs1ErrScaleLength
                        lat1 = latobs(j1)
                        lon1 = lonobs(j1)

        if (i == j1) then
          dist = 0.0
        else
          dist = great_circle_distance(lat2_1, lon2_1, lat1, lon1)
        end if

        if (dist*0.001 < max_dist) then
          b = backErrCov(sigmaBSqr, dist, bErrScaleLength)
          weight = b
          if (i == j1) then
            weight = weight + sigmasqr_obs1
          else
            o1 = obsErrCov(sigmaOSqrj1, oErrScaleLength, dist)
            weight = weight + o1
          end if

          weight = weight * invDataDensities(j1)
          pnew_est(i) = pnew_est(i) + weight * (precep_obs1(j1) - pprev_est(j1))

          weight = b * invDataDensities(j1)
          pnew_ana(i) = pnew_ana(i) + weight * (precep_obs1(j1) - pprev_est(j1))
        end if
      end if

      ! --- Obs2 processing ---
      if (precep_obs2(j1) >= 0.0) then
        sigmaOSqrj1 = sigmasqr_obs2
        oErrScaleLength = obs2ErrScaleLength
        lat1 = latobs(j1) - 0.02d0   ! small offset to avoid collocation
        lon1 = lonobs(j1) - 0.02d0

        if (i == j1) then
          dist = 0.0
        else
          dist = great_circle_distance(lat2_2, lon2_2, lat1, lon1)
        end if

        if (dist*0.001 < max_dist) then
          b = backErrCov(sigmaBSqr, dist, bErrScaleLength)
          weight = b
          if (i == j1) then
            weight = weight + sigmasqr_obs2
          else
            o1 = obsErrCov(sigmaOSqrj1, oErrScaleLength, dist)
            weight = weight + o1
          end if

          weight = weight * invDataDensities(j1)
          pnew_est(i) = pnew_est(i) + weight * (precep_obs2(j1) - pprev_est(j1))

          weight = b * invDataDensities(j1)
          pnew_ana(i) = pnew_ana(i) + weight * (precep_obs2(j1) - pprev_est(j1))
        end if
      end if

    end if
  end do
end do
		end if
      end do

      id = -1   ;       id_incr = 1
      do i=1,npts
         ! Finish analysis and observation estimates for this iteration
            call pick_proc(id,id_incr, numprocs)
            if (id .ne. myid) cycle
            pnew_est(i) = pprev_est(i) + pnew_est(i)
            pnew_ana(i) = pprev_ana(i) + pnew_ana(i)
       end do

      id = -1	;	id_incr = 1
      do i=1,npts
	    call pick_proc(id,id_incr, numprocs)
	    if (id .ne. myid) cycle
            sumObsEstimates(i) = sumObsEstimates(i) + pprev_est(i)
       end do
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
   ! NOTE:  Bratseth values are not interpolated to water points.

   subroutine calc_gridAnalysis(myid,numprocs,npts,obs_typ,sigmasqr_back,sigmasqr_obs1,sigmasqr_obs2,backErrScaleLength,obs1ErrScaleLength,&
                   obs2ErrScaleLength,precep_obs1,precep_obs2,precep_back,&
        sumObsEstimates,npasses,invDataDensities,lat,lon,latobs,lonobs,mrgp_pcp)

      implicit none
      ! Arguments
      ! Arguments
      real*8, allocatable, intent(in) :: lat(:),lon(:),latobs(:),lonobs(:),precep_back(:),precep_obs1(:),precep_obs2(:)
      real*8, allocatable, intent(in) :: sumObsEstimates(:),invDataDensities(:)
      integer, allocatable, intent(in) :: obs_typ(:)
      integer, intent(in) :: npasses,npts,myid, numprocs
      real*8, intent(in) :: sigmasqr_back,sigmasqr_obs1,sigmasqr_obs2,backErrScaleLength,obs1ErrScaleLength,obs2ErrScaleLength
      real*8, allocatable :: mrgp_pcp(:)
      ! Local variables
      real*8 :: tmp_mrgp,weight
      real*8 :: dist,lat1,lon1,sigmaOSqri_1,sigmaOSqri_2,sigmaOSqrj1,sigmaBSqr,oErrScaleLength,bErrScaleLength
      real*8 :: b,max_dist,lat2,lon2
      integer :: i,j,id,x,y,id_incr,i1,j1,xi,yi,xj,yj,jj1,jj2,ii1,ii2,ic
      integer, parameter:: NX=2926, NY=1626

id = -1
id_incr = 1
sigmaBSqr       = sigmasqr_back
bErrScaleLength = backErrScaleLength

! Compute search radius based on available obs scale lengths
if (obs1ErrScaleLength == 0.0) then
    max_dist = 2.0 * MAX(backErrScaleLength, obs2ErrScaleLength) / 1000.0
else if (obs2ErrScaleLength == 0.0) then
    max_dist = 2.0 * MAX(backErrScaleLength, obs1ErrScaleLength) / 1000.0
else
    max_dist = 2.0 * MAX(backErrScaleLength, obs1ErrScaleLength, obs2ErrScaleLength) / 1000.0
end if

ic = int(max_dist / 4.0)

do i = 1, npts
if (precep_back(i) < 0.0) cycle  ! skip invalid points
    call pick_proc(id, id_incr, numprocs)
    if (id .ne. myid) cycle

    tmp_mrgp = precep_back(i)
    if(precep_back(i)>=0.0)  then
    lat2=lat(i)     ;       lon2=lon(i)
        xi  = int((lon(i) - (-169.0)) / 0.04) + 1
        yi  = int((lat(i) -   7.0)   / 0.04) + 1
        ii1 = MAX(1, xi - ic)
        ii2 = MIN(nx, xi + ic)
        jj1 = MAX(1, yi - ic)
        jj2 = MIN(ny, yi + ic)

        do x = ii1, ii2
            do y = jj1, jj2
                j1 = ((x - 1) * ny) + y

                if (precep_back(j1) >= 0.0) then

                    ! ---- Obs1 contribution ----
                    if (precep_obs1(j1) >= 0.0) then
                        lat1 = latobs(j1)
                        lon1 = lonobs(j1)
                        if (i == j1) then
                                dist = 0.0
                        else
                                dist = great_circle_distance(lat2, lon2, lat1, lon1)
                        end if
                        if (dist * 0.001 < max_dist) then
                            b = backErrCov(sigmaBSqr, dist, bErrScaleLength)
                            weight = b * invDataDensities(j1)
                            tmp_mrgp = tmp_mrgp + weight * (npasses * precep_obs1(j1) - sumObsEstimates(j1))
                        end if
                    end if

                    ! ---- Obs2 contribution ----
                    if (precep_obs2(j1) >= 0.0) then
                        lat1 = latobs(j1) - 0.02d0
                        lon1 = lonobs(j1) - 0.02d0
                                if (i == j1) then
                                        dist = 0.0
                                else
                                        dist = great_circle_distance(lat2, lon2, lat1, lon1)
                                end if
                        if (dist * 0.001 < max_dist) then
                            b = backErrCov(sigmaBSqr, dist, bErrScaleLength)
                            weight = b * invDataDensities(j1)
                            tmp_mrgp = tmp_mrgp + weight * (npasses * precep_obs2(j1) - sumObsEstimates(j1))
                        end if
                    end if

                end if
            end do
        end do
 !       write(*,*) i,tmp_mrgp,precep_back(i)
        if (tmp_mrgp < 0.0) tmp_mrgp = 0.0

    end if

    mrgp_pcp(i) = tmp_mrgp

end do
   end subroutine calc_gridAnalysis

  subroutine pick_proc(id, id_incr, numprocs)
    implicit none
    integer,intent(inout) :: id
    integer,intent(inout) :: id_incr
    integer,intent(in) :: numprocs
    id = id + id_incr
    if (id .ge. numprocs) then
       id = numprocs - 1
       id_incr = -1
    else if (id .lt. 0) then
       id = 0
       id_incr = 1
    end if
  end subroutine pick_proc

subroutine read_netcdf(basepath, filename, varname, data, ncid, dimIDs, numLons, numLats)
    character(len=*), intent(in) :: basepath, filename, varname
    integer, intent(out) :: ncid, dimIDs(:), numLons, numLats
    real*8, dimension(:,:), intent(out) :: data
    call check(nf90_open(trim(basepath)//trim(filename), NF90_NOWRITE, ncid))
    call check(nf90_inquire(ncid, ndimsp, nvarsp, nattsp, unlimdimidp))
    call check(nf90_inq_varid(ncid, varname, totprecip_id))
    call check(nf90_inquire_variable(ncid, totprecip_id, dimids=dimIDs))
    call check(nf90_inquire_dimension(ncid, dimIDs(1), len=numLons))
    call check(nf90_inquire_dimension(ncid, dimIDs(2), len=numLats))
    call check(nf90_get_var(ncid, totprecip_id, data))
    call check(nf90_close(ncid))
end subroutine read_netcdf

!-----------------------------------------------------------------------
subroutine create_and_write_netcdf(filename, lat, lon, precep, apcpback, latitude, longitude)
    character(len=*), intent(in) :: filename
    integer, parameter:: NX=2926, NY=1626
    real*8, dimension(NX,NY), intent(in) :: precep, apcpback
    real*8, dimension(NY), intent(in) :: lat
    real*8, dimension(NX), intent(in) :: lon
    real*8, dimension(NX,NY), intent(in) :: latitude, longitude

    integer :: ncid, lat_dimid, lon_dimid
    integer :: lat_id, lon_id, latitude_id, longitude_id
    integer :: totprecip_id, obs_id, apcpobs_id, apcpback_id
    integer, dimension(2) :: dimids1
    type var_info
        character(len=30) :: name
        integer :: id
        character(len=20) :: units
        real :: fill_value
        character(len=50) :: long_name
    end type var_info
    type(var_info), dimension(2) :: vars
    integer :: i

    call check(NF90_create(filename, NF90_clobber, ncid))
    call check(nf90_def_dim(ncid, "latitude", NY, lat_dimid))
    call check(nf90_def_dim(ncid, "longitude", NX, lon_dimid))

    call check(nf90_def_var(ncid, "latitude", NF90_REAL, lat_dimid, lat_id))
    call check(nf90_def_var(ncid, "longitude", NF90_REAL, lon_dimid, lon_id))

    dimids1 = (/ lon_dimid, lat_dimid /)
    call check(nf90_def_var(ncid, "lat", NF90_REAL, dimids1, latitude_id))
    call check(nf90_def_var(ncid, "lon", NF90_REAL, dimids1, longitude_id))

    call check(nf90_put_att(ncid, latitude_id, "units", "degree_north"))
    call check(nf90_put_att(ncid, latitude_id, "_FillValue", -9999.00))
    call check(nf90_put_att(ncid, latitude_id, "long_name", "latitude"))

    call check(nf90_put_att(ncid, longitude_id, "units", "degree_east"))
    call check(nf90_put_att(ncid, longitude_id, "_FillValue", -9999.00))
    call check(nf90_put_att(ncid, longitude_id, "long_name", "longitude"))
    ! Define the variables first, get their IDs
    call check(nf90_def_var(ncid, "TotalPrecip_tavg", NF90_REAL, dimids1, totprecip_id))
    call check(nf90_def_var(ncid, "PrecepBack", NF90_REAL, dimids1, apcpback_id))

    ! Assign the vars array with now valid IDs
    vars(1) = var_info("TotalPrecip_tavg", totprecip_id, "mm", -9999.00, "total precipitation")
    vars(2) = var_info("PrecepBack", apcpback_id, "mm", -9999.00, "background precipitation")
    do i = 1, size(vars)
        call check(nf90_put_att(ncid, vars(i)%id, "units", vars(i)%units))
        if (vars(i)%fill_value /= 0.0) then
            call check(nf90_put_att(ncid, vars(i)%id, "_FillValue", vars(i)%fill_value))
        end if
        call check(nf90_put_att(ncid, vars(i)%id, "long_name", vars(i)%long_name))
    end do
    call check(nf90_enddef(ncid))

    call check(nf90_put_var(ncid, lat_id, lat))
    call check(nf90_put_var(ncid, lon_id, lon))
    call check(nf90_put_var(ncid, totprecip_id, precep))
    call check(nf90_put_var(ncid, apcpback_id, apcpback))
    call check(nf90_put_var(ncid, latitude_id, latitude))
    call check(nf90_put_var(ncid, longitude_id, longitude))
    call check(nf90_close(ncid))
end subroutine create_and_write_netcdf

   !---------------------------------------------------------------------------
   ! Observation error covariance function.
   real*8 function obsErrCov(sigmaOSqr,oErrScaleLength,dist)
      implicit none
      real*8, intent(in) :: sigmaOSqr
      real*8, intent(in) :: oErrScaleLength ! in meters
      real*8, intent(in) :: dist ! in meters
      obsErrCov = sigmaOSqr*obsErrCorr(oErrScaleLength,dist)
   end function obsErrCov

   !---------------------------------------------------------------------------
   ! Observation error correlation function.  Currently Gaussian.
   real*8 function obsErrCorr(oErrScaleLength,dist)
      implicit none
      real*8, intent(in) :: oErrScaleLength ! in meters
      real*8, intent(in) :: dist ! in meters
      real*8 :: invOErrScaleLength
      invOErrScaleLength = 1. / oErrScaleLength
      obsErrCorr = exp(-1*dist*dist*invOErrScaleLength*InvOErrScaleLength)
   end function obsErrCorr

   !---------------------------------------------------------------------------
   ! Background error covariance function.
   real*8 function backErrCov(sigmaBSr,dist,scale_length)
      implicit none
      real*8, intent(in) :: sigmaBSr
      real*8, intent(in) :: dist ! in meters
      real*8, intent(in) :: scale_length
      backErrCov = sigmaBSr*backErrCorr(dist,scale_length)
   end function backErrCov

   !---------------------------------------------------------------------------
   ! Background error correlation function.  Currently Gaussian.
   real*8 function backErrCorr(dist,scale_length)
      implicit none
      real*8, intent(in) :: dist ! in meters
      real*8, intent(in) :: scale_length
      real*8 :: inv_scale_length
      inv_scale_length = 1./scale_length
      backErrCorr = exp(-1*dist*dist*inv_scale_length*inv_scale_length)
   end function backErrCorr

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
