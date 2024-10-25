!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------

module USAF_readFilesMod

  implicit none
  private

  public :: readfiles

contains

  subroutine readfiles(myid,numprocs, starttime, endtime, deltatime, &
        datadir, max_stations, max_sat_reports, max_vario_bins, &
        vario_bin_dist, vario, icounts_vario, &
        imax_lon, jmax_lat, dlon, dlat, dist_thresh, &
        is_sattype, &
        use_blacklist, nstns, blacklist_stns)

    ! Imports
    use esmf
    use mpi
    use USAF_GridHashMod, only: GridHash, newGridHash, destroyGridHash, &
         insertIntoGridHash, getObindexVectorFromGridHash, &
         createIJForGridHash
    use USAF_ReportsMod, only: Reports, newReports, getNobs, getReport, &
         destroyReports, appendToReports, bcast_reports
    use USAF_StationsMod, only: great_circle_distance
    use USAF_sharedMod, only: is_gage, endrun, pick_proc
    use USAF_StnOBDictMod, only: stnOBDict_t

    ! Defaults
    implicit none

    ! Arguments
    integer, intent(in) :: myid
    integer, intent(in) :: numprocs
    type(esmf_time), intent(in) :: starttime
    type(esmf_time), intent(in) :: endtime
    type(esmf_timeinterval), intent(in) :: deltatime
    character(len=255), intent(in) :: datadir
    integer, intent(in) :: max_stations
    integer, intent(in) :: max_sat_reports
    integer, intent(in) :: max_vario_bins
    integer, intent(in) :: vario_bin_dist
    double precision, intent(inout) :: vario(max_vario_bins)
    integer*8, intent(inout) :: icounts_vario(max_vario_bins)
    integer, intent(in) :: imax_lon
    integer, intent(in) :: jmax_lat
    real, intent(in) :: dlon
    real, intent(in) :: dlat
    real, intent(in) :: dist_thresh
    logical, external :: is_sattype ! Function
    logical, intent(in) :: use_blacklist
    integer, intent(in) :: nstns
    character(len=9), allocatable, intent(in) ::  blacklist_stns(:)

    ! Locals
    type(StnOBDIct_t) :: obStnHash
    integer :: iunit_input
    integer*8 :: sampleSize
    type(esmf_time) :: curtime
    integer :: yyyy, mm, dd, hh
    integer :: rc
    character(len=44) :: infile
    character(len=123) :: fullpath
    integer :: istat
    integer :: ierr
    integer :: count_skips
    type(Reports) :: R_gages_sat, R_matches
    integer :: icount_gages, icount_sat
    character(len=80) :: line
    character(len=10) :: network, network_k, network_kk, network_new
    character(len=10) :: platform, platform_k, platform_kk
    real :: latitude, longitude, lat_k, lon_k, lat_kk, lon_kk
    real :: O, B, A, O_k, B_k, O_kk, B_kk, B_new, OMB_k, OMB_kk
    logical :: skip
    integer :: istn
    character(len=10) :: yyyymmddhh, yyyymmddhh_k, yyyymmddhh_kk
    integer :: nobs
    type(GridHash) :: sat2d, gages2d
    integer :: i, j, ii, jj, k, kk
    integer :: nobs_gages
    integer, allocatable :: gage_obindexVector(:)
    integer, allocatable :: sat_obindexVector(:)
    real :: min_dist_k_kk
    integer :: jdelta_lat, idelta_lon
    integer :: nobs_sat
    real :: dist_k_kk
    character(len=10), allocatable :: stationNames(:)
    real, allocatable :: meanOMBS(:)
    integer :: numKeys
    double precision, allocatable :: vario_proc(:)
    integer*4, allocatable :: icounts_vario_proc_i4(:)
    double precision, allocatable :: vario_allproc(:)
    integer*4, allocatable :: icounts_vario_allproc_i4(:)
    integer :: sample_size_proc_i4, sample_size_allproc_i4
    integer :: id, id_incr
    integer :: index
    character(len=10) :: c_k

    iunit_input = 10
    sampleSize = 0

    !call ESMF_LogWrite( &
    !     "Creating obStnHash", &
    !     ESMF_LOGMSG_INFO)
    !call ESMF_LogFlush()

    call obStnHash%new()

    !call ESMF_LogWrite( &
    !     "Creating R_matches", &
    !     ESMF_LOGMSG_INFO)
    !call ESMF_LogFlush()

    R_matches = newReports(28*4*max_stations)

    ! We will loop through each file, saving data as appropriate
    curtime = starttime
    do
       if (curtime .ge. endtime) exit
       call esmf_timeget(curtime, &
            yy=yyyy, mm=mm, dd=dd, h=hh, rc=rc)
       if (rc .ne. ESMF_SUCCESS) then
          call ESMF_LogWrite( &
               "Cannot get current time!", &
               ESMF_LOGMSG_ERROR)
          call endrun(1)
       end if

       ! Advance curtime to next time so we can cycle if a problem occurs
       curtime = curtime + deltatime

       if (myid .eq. 0) then
          write(yyyymmddhh, '(I4.4,I2.2,I2.2,I2.2)') yyyy, mm, dd, hh
          write(infile, 1000) 'oba_', yyyy, mm, dd, hh, '_12.txt'
1000      format(A,I4.4,I2.2,I2.2,I2.2,A)
          fullpath = trim(datadir) // '/' // trim(infile)
          open(unit=iunit_input, &
               file=trim(fullpath), &
               status="OLD", &
               iostat=istat)
       end if
       call mpi_bcast(istat, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
       if (ierr .ne. MPI_SUCCESS) then
          call ESMF_LogWrite( &
               "Problem calling mpi_bcast!", &
               ESMF_LOGMSG_ERROR)
          call endrun(1)
       end if
       if (istat .ne. 0) then
          if (myid .eq. 0) then
             call ESMF_LogWrite( &
                  'Problem opening '//trim(fullpath), &
                  ESMF_LOGMSG_WARNING)
             close(unit=iunit_input)
          end if
          cycle
       end if

       if (myid .eq. 0) then
          call ESMF_LogWrite( &
               'Reading file '//trim(infile), &
               ESMF_LOGMSG_INFO)
          call ESMF_LogFlush()
       end if

       ! First, get all the gages and specified satellite obs in the
       ! current file.
       count_skips = 0

       ! call ESMF_LogWrite( &
       !      "Creating R_gages_sat", &
       !      ESMF_LOGMSG_INFO)
       ! call ESMF_LogFlush()

       R_gages_sat = newReports(max_stations+max_sat_reports)
       if (myid .eq. 0) then
          icount_gages = 0
          icount_sat = 0
          do
             read(iunit_input, '(A)', iostat=istat) line
             if (istat .ne. 0) exit
             if (line(2:2) .eq. '#') cycle
             read(line, &
                  '(A10,1x,A10,1x,F8.3,1x,F8.3,1x,F8.3,1x,F8.3,1x,F8.3)') &
                  network, platform, latitude, longitude, O, B, A
             if (is_gage(adjustl(network))) then

                ! Apply blacklist
                skip = .false.
                if (use_blacklist) then
                   do istn = 1, nstns
                      if (trim(blacklist_stns(istn)) .eq. &
                           adjustl(platform)) then
                         count_skips = count_skips + 1
                         skip = .true.
                         exit
                      end if
                   end do
                end if
                if (skip) cycle

                icount_gages = icount_gages + 1
                call appendToReports(R_gages_sat, network, platform, &
                     yyyymmddhh, &
                     latitude, longitude, O, B, A)
             else if (is_sattype(adjustl(network))) then
                icount_sat = icount_sat + 1
                call appendToReports(R_gages_sat, network, platform,&
                     yyyymmddhh, &
                     latitude, longitude, O, B, A)
             else
                cycle
             end if
          end do
          close(unit=iunit_input)
       end if

       ! Share reports across processors
       call mpi_barrier(mpi_comm_world, ierr)
       call bcast_reports(R_gages_sat, myid, ierr)
       if (ierr .ne. MPI_SUCCESS) then
          call ESMF_LogWrite( &
               'Cannot share reports across processors!', &
               ESMF_LOGMSG_ERROR)
          call endrun(1)
       end if

       ! Create geographic lookup tables for all data.  We do this here
       ! since the hash tables use linked lists under the hood, which
       ! are difficult to pass via MPI.
       call mpi_barrier(mpi_comm_world, ierr)
       nobs = getNobs(R_gages_sat)

!       call ESMF_LogWrite( &
!            "Creating sat2d", &
!            ESMF_LOGMSG_INFO)
!       call ESMF_LogFlush()

       sat2d = newGridHash(IMAX_LON, JMAX_LAT, DLON, DLAT)

!       call ESMF_LogWrite( &
!            "Creating gages2d", &
!            ESMF_LOGMSG_INFO)
!       call ESMF_LogFlush()

       gages2d = newGridHash(IMAX_LON, JMAX_LAT, DLON, DLAT)
       do k = 1, nobs
          call getReport(R_gages_sat, k, network=network_k, &
               latitude=lat_k, longitude=lon_k)
          if (is_gage(adjustl(network_k))) then
             call createIJForGridHash(gages2d, lat_k, lon_k, ii, jj)
             call insertIntoGridHash(gages2d, ii, jj, k)
          else
             call createIJForGridHash(sat2d, lat_k, lon_k, ii, jj)
             call insertIntoGridHash(sat2d, ii, jj, k)
          end if
       end do ! k

!       call ESMF_LogWrite( &
!            "HERE", &
!            ESMF_LOGMSG_INFO)
!       call ESMF_LogFlush()

       ! Now we try to match each gage to the nearest sat ob
       call mpi_barrier(mpi_comm_world, ierr)
       do j = 1, JMAX_LAT
          do i = 1, IMAX_LON
             ! Get list of gages in this box
             call getObindexVectorFromGridHash(gages2d, i, j, &
                  nobs_gages, gage_obindexVector)
             if (nobs_gages .eq. 0) cycle

             do k = 1, nobs_gages

                call getReport(R_gages_sat, gage_obindexVector(k), &
                     network=network_k, &
                     platform=platform_k, yyyymmddhh=yyyymmddhh_k,&
                     latitude=lat_k, longitude=lon_k, O=O_k, B=B_k)
                min_dist_k_kk = 9999999999999999.
                ! Now search 3x3 boxes for nearest sat ob
                do jdelta_lat = -1, 1
                   jj = j + jdelta_lat
                   if (jj .lt. 1) cycle
                   if (jj .gt. JMAX_LAT) cycle
                   do idelta_lon = -1, 1
                      ii = i + idelta_lon
                      if (ii .lt. 1) then
                         ii = ii + IMAX_LON
                      else if (ii .gt. IMAX_LON) then
                         ii = ii - IMAX_LON
                      end if

                      ! Get list of sat report in this box
                      call getObindexVectorFromGridHash(sat2d, ii, jj, &
                           nobs_sat, sat_obindexVector)
                      if (nobs_sat .eq. 0) cycle

                      do kk = 1, nobs_sat
                         call getReport(R_gages_sat, &
                              sat_obindexVector(kk), network=network_kk, &
                              platform=platform_kk, &
                              yyyymmddhh=yyyymmddhh_kk,&
                              latitude=lat_kk, longitude=lon_kk, O=O_kk)
                         dist_k_kk = &
                              great_circle_distance(lat_k, lon_k, &
                              lat_kk, lon_kk)
                         if (dist_k_kk .lt. min_dist_k_kk) then
                            min_dist_k_kk = dist_k_kk
                            B_new = O_kk ! Use sat retrieval as background
                         end if
                      end do ! kk
                      if (allocated(sat_obindexVector)) then
                         deallocate(sat_obindexVector)
                      end if
                   end do ! idelta_lon
                end do ! idelta_lat

                ! See if nearest sat is close enough to gage
                if (min_dist_k_kk .lt. dist_thresh) then
                   network_new = "SATGAGE"
                   call appendToReports(R_matches,&
                        network=network_new, &
                        platform=platform_k, &
                        yyyymmddhh=yyyymmddhh_k,&
                        latitude=lat_k, &
                        longitude=lon_k, &
                        O=O_k,&
                        B=B_new,&
                        A=0.)
!                   call ESMF_LogWrite( &
!                        "Adding platform_k: "//trim(platform_k), &
!                        ESMF_LOGMSG_INFO)
!                   call ESMF_LogFlush()

                   call obStnHash%insert(platform_k, O_k, B_new)
                end if
             end do ! k
             if (allocated(gage_obindexVector)) then
                deallocate(gage_obindexVector)
             end if
          end do ! ilon
       end do ! jlat

!       call ESMF_LogWrite( &
!            "destroying R_gages_sat", &
!            ESMF_LOGMSG_INFO)
!       call ESMF_LogFlush()

       ! Clean up
       call destroyReports(R_gages_sat)

!       call ESMF_LogWrite( &
!            "destroying gages2d", &
!            ESMF_LOGMSG_INFO)
!       call ESMF_LogFlush()

       call destroyGridHash(gages2d)

!       call ESMF_LogWrite( &
!            "destroying sat2d", &
!            ESMF_LOGMSG_INFO)
!       call ESMF_LogFlush()

       call destroyGridHash(sat2d)

    enddo ! next file

    ! We have all the matches from the files.  Now, calculate the mean OMBs for
    ! each gage-satellite pair.

    call ESMF_LogWrite( &
         "Calling getMeanOMBs", &
         ESMF_LOGMSG_INFO)
    call ESMF_LogFlush()
    call obStnHash%getMeanOMBs(numKeys, stationNames, meanOMBs)

    call ESMF_LogWrite( &
         "Back from getMeanOMBs", &
         ESMF_LOGMSG_INFO)
    call ESMF_LogFlush()


    ! Now, update semivariogram sums and counts, excluding those stations with
    ! large meanOMBs
    nobs = getNobs(R_matches)

    write(c_k, '(I10)') nobs
    call ESMF_LogWrite("EMK: k = "//c_k, &
         ESMF_LOGMSG_INFO)
    call ESMF_LogFlush()

    allocate(vario_proc(max_vario_bins))
    allocate(icounts_vario_proc_i4(max_vario_bins))
    do k = 1, nobs
       vario_proc = 0
       icounts_vario_proc_i4 = 0
       id = -1
       id_incr = 1
       sample_size_proc_i4 = 0
       call getReport(R_matches, k, platform=platform_k,&
            yyyymmddhh=yyyymmddhh_k, &
            latitude=lat_k, longitude=lon_k, &
            O=O_k, B=B_k)
       skip = .false.
       do i = 1, numKeys
          if (trim(platform_k) .eq. trim(stationNames(i))) then
             if (meanOMBs(i) > 10) then
                skip = .true.
                exit
             end if
          end if
       end do
       if (skip) cycle

       write(c_k, '(I5)') k
       call ESMF_LogWrite("EMK: k = "//c_k, &
            ESMF_LOGMSG_INFO)
       call ESMF_LogFlush()

       OMB_k = O_k - B_k

       do kk = k+1, nobs

          ! See if current processor is responsible for this j ob.
          call pick_proc(id, id_incr, numprocs)
          if (id .ne. myid) cycle

          call getReport(R_matches, kk, platform=platform_kk,&
               yyyymmddhh=yyyymmddhh_kk, &
               latitude=lat_kk, longitude=lon_kk, &
               O=O_kk, B=B_kk)
          if (trim(yyyymmddhh_k) .ne. trim(yyyymmddhh_kk)) then
             exit ! Subsequent kk reports will be a different time.
          end if
          skip = .false.
          do i = 1, numKeys
             if (trim(platform_k) .eq. trim(stationNames(i))) then
                if (meanOMBs(i) > 10) then
                   skip = .true.
                   exit
                end if
             end if
          end do ! i
          if (skip) cycle

          OMB_kk = O_kk - B_kk

          dist_k_kk = great_circle_distance(lat_k, lon_k, lat_kk, lon_kk)
          index = int(dist_k_kk*0.001/vario_bin_dist) + 1
          if (index .gt. max_vario_bins) cycle

          vario_proc(index) = vario_proc(index) + &
               ((OMB_k-OMB_kk)*(OMB_k-OMB_kk))
          icounts_vario_proc_i4(index) = icounts_vario_proc_i4(index) + 1
          sample_size_proc_i4 = sample_size_proc_i4 + 1
       end do ! kk

       ! Now need to collect the data from each processor
       sample_size_allproc_i4 = 0
       call mpi_barrier(mpi_comm_world, ierr)
       call mpi_allreduce(sample_size_proc_i4, sample_size_allproc_i4, &
            1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
       sampleSize = sampleSize + sample_size_allproc_i4

       allocate(vario_allproc(max_vario_bins))
       vario_allproc = 0
       call mpi_barrier(mpi_comm_world,ierr)
       call mpi_allreduce(vario_proc, vario_allproc, max_vario_bins, &
            MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
       vario = vario + vario_allproc

       allocate(icounts_vario_allproc_i4(max_vario_bins))
       icounts_vario_allproc_i4 = 0
       call mpi_barrier(mpi_comm_world, ierr)
       call mpi_allreduce(icounts_vario_proc_i4, &
            icounts_vario_allproc_i4, &
            max_vario_bins, &
            MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
       icounts_vario = icounts_vario + icounts_vario_allproc_i4
       deallocate(icounts_vario_allproc_i4)
       deallocate(vario_allproc)
    end do ! k

    deallocate(vario_proc)
    deallocate(icounts_vario_proc_i4)


    call destroyReports(R_matches)
    call obStnHash%destroy()
    deallocate(stationNames)
    deallocate(meanOMBs)

  end subroutine readfiles

end module USAF_readFilesMod
