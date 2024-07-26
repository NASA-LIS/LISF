!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------

!
! ROUTINE: procOBA_Sat
!
! DESCRIPTION:
! Main program for generating semivariogram of gage-satellite differences
! with distance.
!
! REVISION HISTORY:
! 26 0ct 2020:  Eric Kemp.  Initial Specification.  Multi-process
!   MPI doesn't work yet, just use single process.
!
program main

   ! Imports
   use esmf
   use mpi
   use USAF_sharedMod

   ! Defaults
   implicit none

   ! Local variables
   character(len=255) :: outfile
   integer :: iunit_input, istat, iunit_out_vario
   integer :: j, i
   integer*8, allocatable :: icounts_vario(:)
   double precision, allocatable :: vario(:)
   character(ESMF_MAXPATHLEN) :: cfgfile
   type(ESMF_Config) :: cf
   integer :: rc
   integer :: max_stations
   integer :: max_sat_reports
   integer :: max_num_files
   integer :: max_vario_bins
   integer :: vario_bin_dist
   character(len=255) :: datadir
   integer :: startyear
   integer :: startmonth
   integer :: startday
   integer :: starthour
   integer :: endyear
   integer :: endmonth
   integer :: endday
   integer :: endhour
   integer :: intervalyear
   integer :: intervalmonth
   integer :: intervalday
   integer :: intervalhour
   type(ESMF_Time) :: starttime, endtime, curtime
   type(ESMF_TimeInterval) :: deltatime
   type(ESMF_VM) :: vm
   integer :: icount
   integer, parameter :: maxlen_sattype = 9
   character(len=maxlen_sattype) :: sattype
   integer, parameter :: max_sattypes = 4
   character(len=maxlen_sattype) :: sattypes(max_sattypes)
   logical :: found_sattype
   integer :: numprocs, ierr
   integer :: myid
   integer :: num_args
   integer :: imax_lon
   integer :: jmax_lat
   real :: dlon
   real :: dlat
   real :: dist_thresh
   logical :: use_blacklist
   character(len=255) :: blacklist_file
   character(len=9), allocatable :: blacklist_stns(:)
   integer :: nstns
   character(len=255) :: logname

   character(len=maxlen_sattype), parameter :: SSMI      = "SSMI"
   character(len=maxlen_sattype), parameter :: CMORPH    = "CMORPH"
   character(len=maxlen_sattype), parameter :: GEOPRECIP = "GEOPRECIP"
   character(len=maxlen_sattype), parameter :: IMERG     = "IMERG"
   sattypes = (/SSMI, CMORPH, GEOPRECIP, IMERG/)

   ! Set logfile name
   logname = 'procOBA_Sat.log'
   num_args = command_argument_count()
   if (num_args .eq. 2) then
      call get_command_argument(2, logname)
   end if

   ! Initialize ESMF.  Must happen first.  This calls MPI_Init under the hood.
   call esmf_initialize(vm=vm, defaultCalKind=ESMF_CALKIND_GREGORIAN, &
        defaultLogFileName=trim(logname), rc=rc)
   if (rc .ne. ESMF_SUCCESS) then
      call ESMF_LogWrite("Cannot initialize ESMF!", ESMF_LOGMSG_ERROR)
      call endrun(1)
   end if

   ! Get some MPI information
   call mpi_comm_rank(MPI_COMM_WORLD, myid, ierr)
   if (ierr .ne. MPI_SUCCESS) then
      call ESMF_LogWrite("Problem calling mpi_comm_rank", &
           ESMF_LOGMSG_ERROR)
      call endrun(1)
   end if
   call mpi_comm_size(MPI_COMM_WORLD, numprocs, ierr)
   if (ierr .ne. MPI_SUCCESS) then
      call ESMF_LogWrite("Problem calling mpi_comm_size", &
           ESMF_LOGMSG_ERROR)
      call endrun(1)
   end if

   if (myid .eq. 0) then
      num_args = command_argument_count()
      if (num_args .ne. 2) then
         call ESMF_LogWrite("Improper program invocation", &
              ESMF_LOGMSG_ERROR)
         call ESMF_LogWrite("USAGE: procOBA_NWP <configfile> <logfile>", &
              ESMF_LOGMSG_ERROR)
         call endrun(1)
      end if
      call get_command_argument(1, cfgfile)
   end if

   if (myid .eq. 0) then
      ! Read config file
      cf = ESMF_ConfigCreate(rc=rc)
      call esmf_ConfigLoadFile(cf, cfgfile, rc=rc)
      if (rc .ne. ESMF_SUCCESS) then
         call ESMF_LogWrite("Cannot open "//trim(cfgfile), &
              ESMF_LOGMSG_ERROR)
         call endrun(1)
      end if
   end if

   max_stations = 0
   if (myid .eq. 0) then
      call esmf_configgetattribute(cf, max_stations, label="max_stations:", &
           rc=rc)
      if (rc .ne. ESMF_SUCCESS .or. max_stations .lt. 1) then
         call ESMF_LogWrite( &
              "Must specify positive max_stations in config file!", &
              ESMF_LOGMSG_ERROR)
         call endrun(1)
      end if
   end if
   call mpi_barrier(mpi_comm_world, ierr)
   call mpi_bcast(max_stations, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
   if (ierr .ne. MPI_SUCCESS) then
      call ESMF_LogWrite( &
           "Problem calling mpi_bcast!", &
           ESMF_LOGMSG_ERROR)
      call endrun(1)
   end if

   max_sat_reports = 0
   if (myid .eq. 0) then
      call esmf_configgetattribute(cf, max_sat_reports,&
           label="max_sat_reports:", rc=rc)
      if (rc .ne. ESMF_SUCCESS .or. max_sat_reports .lt. 1) then
         call ESMF_LogWrite( &
              "Must specify max_sat_reports in config file!", &
              ESMF_LOGMSG_ERROR)
         call endrun(1)
      end if
   end if
   call mpi_barrier(mpi_comm_world, ierr)
   call mpi_bcast(max_sat_reports, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
   if (ierr .ne. MPI_SUCCESS) then
      call ESMF_LogWrite( &
           "Problem calling mpi_bcast!", &
           ESMF_LOGMSG_ERROR)
      call endrun(1)
   end if

   max_vario_bins = 0
   if (myid .eq. 0) then
      call esmf_configgetattribute(cf, max_vario_bins, &
           label="max_vario_bins:", &
           rc=rc)
      if (rc .ne. ESMF_SUCCESS .or. max_vario_bins .lt. 1) then
         call ESMF_LogWrite( &
              "Must specify positive max_vario_bins in config file!", &
              ESMF_LOGMSG_ERROR)
         call endrun(1)
      end if
   end if
   call mpi_barrier(mpi_comm_world, ierr)
   call mpi_bcast(max_vario_bins, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
   if (ierr .ne. MPI_SUCCESS) then
      call ESMF_LogWrite( &
           "Problem calling mpi_bcast!", &
           ESMF_LOGMSG_ERROR)
      call endrun(1)
   end if

   vario_bin_dist = 0
   if (myid .eq. 0) then
      call esmf_configgetattribute(cf, vario_bin_dist, &
           label="vario_bin_dist:", &
           rc=rc)
      if (rc .ne. ESMF_SUCCESS .or. vario_bin_dist .lt. 1) then
         call ESMF_LogWrite( &
              "Must specify positive vario_bin_dist in config file!", &
              ESMF_LOGMSG_ERROR)
         call endrun(1)
      end if
   end if
   call mpi_barrier(mpi_comm_world, ierr)
   call mpi_bcast(vario_bin_dist, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
   if (ierr .ne. MPI_SUCCESS) then
      call ESMF_LogWrite( &
           "Problem calling mpi_bcast!", &
           ESMF_LOGMSG_ERROR)
      call endrun(1)
   end if

   datadir = 'NULL'
   if (myid .eq. 0) then
      call esmf_configgetattribute(cf, datadir, &
           label="datadir:", &
           rc=rc)
      if (rc .ne. ESMF_SUCCESS) then
         call ESMF_LogWrite( &
              "Must specify datadir in config file!", &
              ESMF_LOGMSG_ERROR)
         call endrun(1)
      end if
   end if
   call mpi_barrier(mpi_comm_world, ierr)
   call mpi_bcast(datadir, len(datadir), MPI_CHARACTER, 0, MPI_COMM_WORLD, &
        ierr)
   if (ierr .ne. MPI_SUCCESS) then
      call ESMF_LogWrite( &
           "Problem calling mpi_bcast!", &
           ESMF_LOGMSG_ERROR)
      call endrun(1)
   end if

   startyear = 0
   if (myid .eq. 0) then
      call esmf_configgetattribute(cf, startyear, &
           label="startyear:", &
           rc=rc)
      if (rc .ne. ESMF_SUCCESS) then
         call ESMF_LogWrite( &
              "Must specify startyear in config file!", &
              ESMF_LOGMSG_ERROR)
         call endrun(1)
      end if
   end if
   call mpi_barrier(mpi_comm_world, ierr)
   call mpi_bcast(startyear, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
   if (ierr .ne. MPI_SUCCESS) then
      call ESMF_LogWrite( &
           "Problem calling mpi_bcast!", &
           ESMF_LOGMSG_ERROR)
      call endrun(1)
   end if

   startmonth = 0
   if (myid .eq. 0) then
      call esmf_configgetattribute(cf, startmonth, &
           label="startmonth:", &
           rc=rc)
      if (rc .ne. ESMF_SUCCESS) then
         call ESMF_LogWrite( &
              "Must specify startmonth in config file!", &
              ESMF_LOGMSG_ERROR)
         call endrun(1)
      end if
   end if
   call mpi_barrier(mpi_comm_world, ierr)
   call mpi_bcast(startmonth, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
   if (ierr .ne. MPI_SUCCESS) then
         call ESMF_LogWrite( &
              "Problem calling mpi_bcast!", &
              ESMF_LOGMSG_ERROR)
      call endrun(1)
   end if

   startday = 0
   if (myid .eq. 0) then
      call esmf_configgetattribute(cf, startday, &
           label="startday:", &
           rc=rc)
      if (rc .ne. ESMF_SUCCESS) then
         call ESMF_LogWrite( &
              "Must specify startday in config file!", &
              ESMF_LOGMSG_ERROR)
         call endrun(1)
      end if
   end if
   call mpi_barrier(mpi_comm_world, ierr)
   call mpi_bcast(startday, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
   if (ierr .ne. MPI_SUCCESS) then
      call ESMF_LogWrite( &
           "Problem calling mpi_bcast!", &
           ESMF_LOGMSG_ERROR)
      call endrun(1)
   end if

   starthour = 0
   if (myid .eq. 0) then
      call esmf_configgetattribute(cf, starthour, &
           label="starthour:", &
           rc=rc)
      if (rc .ne. ESMF_SUCCESS) then
         call ESMF_LogWrite( &
              "Must specify starthour in config file!", &
              ESMF_LOGMSG_ERROR)
         call endrun(1)
      end if
   end if
   call mpi_barrier(mpi_comm_world, ierr)
   call mpi_bcast(starthour, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
   if (ierr .ne. MPI_SUCCESS) then
      call ESMF_LogWrite( &
           "Problem calling mpi_bcast!", &
           ESMF_LOGMSG_ERROR)
      call endrun(1)
   end if

   endyear = 0
   if (myid .eq. 0) then
      call esmf_configgetattribute(cf, endyear, &
           label="endyear:", &
           rc=rc)
      if (rc .ne. ESMF_SUCCESS) then
         call ESMF_LogWrite( &
              "Must specify endyear in config file!", &
              ESMF_LOGMSG_ERROR)
         call endrun(1)
      end if
   end if
   call mpi_barrier(mpi_comm_world, ierr)
   call mpi_bcast(endyear, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
   if (ierr .ne. MPI_SUCCESS) then
      call ESMF_LogWrite( &
           "Problem calling mpi_bcast!", &
           ESMF_LOGMSG_ERROR)
      call endrun(1)
   end if

   endmonth = 0
   if (myid .eq. 0) then
      call esmf_configgetattribute(cf, endmonth, &
           label="endmonth:", &
           rc=rc)
      if (rc .ne. ESMF_SUCCESS) then
         call ESMF_LogWrite( &
              "Must specify endmonth in config file!", &
              ESMF_LOGMSG_ERROR)
         call endrun(1)
      end if
   end if
   call mpi_barrier(mpi_comm_world, ierr)
   call mpi_bcast(endmonth, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
   if (ierr .ne. MPI_SUCCESS) then
      call ESMF_LogWrite( &
           "Problem calling mpi_bcast!", &
           ESMF_LOGMSG_ERROR)
      call endrun(1)
   end if

   endday = 0
   if (myid .eq. 0) then
      call esmf_configgetattribute(cf, endday, &
           label="endday:", &
           rc=rc)
      if (rc .ne. ESMF_SUCCESS) then
         call ESMF_LogWrite( &
              "Must specify endday in config file!", &
              ESMF_LOGMSG_ERROR)
         call endrun(1)
      end if
   end if
   call mpi_barrier(mpi_comm_world, ierr)
   call mpi_bcast(endday, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
   if (ierr .ne. MPI_SUCCESS) then
      call ESMF_LogWrite( &
           "Problem calling mpi_bcast!", &
           ESMF_LOGMSG_ERROR)
      call endrun(1)
   end if

   endhour = 0
   if (myid .eq. 0) then
      call esmf_configgetattribute(cf, endhour, &
           label="endhour:", &
           rc=rc)
      if (rc .ne. ESMF_SUCCESS) then
         call ESMF_LogWrite( &
              "Must specify endhour in config file!", &
              ESMF_LOGMSG_ERROR)
         call endrun(1)
      end if
   end if
   call mpi_barrier(mpi_comm_world, ierr)
   call mpi_bcast(endhour, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
   if (ierr .ne. MPI_SUCCESS) then
      call ESMF_LogWrite( &
           "Problem calling mpi_bcast!", &
           ESMF_LOGMSG_ERROR)
      call endrun(1)
   end if

   intervalyear = 0
   if (myid .eq. 0) then
      call esmf_configgetattribute(cf, intervalyear, &
           label="intervalyear:", &
           rc=rc)
      if (rc .ne. ESMF_SUCCESS) then
         call ESMF_LogWrite( &
              "Must specify intervalyear in config file!", &
              ESMF_LOGMSG_ERROR)
         call endrun(1)
      end if
   end if
   call mpi_barrier(mpi_comm_world, ierr)
   call mpi_bcast(intervalyear, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
   if (ierr .ne. MPI_SUCCESS) then
      call ESMF_LogWrite( &
           "Problem calling mpi_bcast!", &
           ESMF_LOGMSG_ERROR)
      call endrun(1)
   end if

   intervalmonth = 0
   if (myid .eq. 0) then
      call esmf_configgetattribute(cf, intervalmonth, &
           label="intervalmonth:", &
           rc=rc)
      if (rc .ne. ESMF_SUCCESS) then
         call ESMF_LogWrite( &
              "Must specify intervalmonth in config file!", &
              ESMF_LOGMSG_ERROR)
         call endrun(1)
      end if
   end if
   call mpi_barrier(mpi_comm_world, ierr)
   call mpi_bcast(intervalmonth, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
   if (ierr .ne. MPI_SUCCESS) then
      call ESMF_LogWrite( &
           "Problem calling mpi_bcast!", &
           ESMF_LOGMSG_ERROR)
      call endrun(1)
   end if

   intervalday = 0
   if (myid .eq. 0) then
      call esmf_configgetattribute(cf, intervalday, &
           label="intervalday:", &
           rc=rc)
      if (rc .ne. ESMF_SUCCESS) then
         call ESMF_LogWrite( &
              "Must specify intervalday in config file!", &
              ESMF_LOGMSG_ERROR)
         call endrun(1)
      end if
   end if
   call mpi_barrier(mpi_comm_world, ierr)
   call mpi_bcast(intervalday, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
   if (ierr .ne. MPI_SUCCESS) then
      call ESMF_LogWrite( &
           "Problem calling mpi_bcast!", &
           ESMF_LOGMSG_ERROR)
      call endrun(1)
   end if

   intervalhour = 0
   if (myid .eq. 0) then
      call esmf_configgetattribute(cf, intervalhour, &
           label="intervalhour:", &
           rc=rc)
      if (rc .ne. ESMF_SUCCESS) then
         call ESMF_LogWrite( &
              "Must specify intervalhour in config file!", &
              ESMF_LOGMSG_ERROR)
         call endrun(1)
      end if
   end if
   call mpi_barrier(mpi_comm_world, ierr)
   call mpi_bcast(intervalhour, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
   if (ierr .ne. MPI_SUCCESS) then
      call ESMF_LogWrite( &
           "Problem calling mpi_bcast!", &
           ESMF_LOGMSG_ERROR)
      call endrun(1)
   end if

   outfile = "NULL"
   if (myid .eq. 0) then
      call esmf_configgetattribute(cf, outfile, &
           label="outfile:", &
           rc=rc)
      if (rc .ne. ESMF_SUCCESS) then
         call ESMF_LogWrite( &
              "Must specify outfile in config file!", &
              ESMF_LOGMSG_ERROR)
         call endrun(1)
      end if
   end if
   call mpi_barrier(mpi_comm_world, ierr)
   call mpi_bcast(outfile, len(outfile), MPI_CHARACTER, 0, MPI_COMM_WORLD, &
        ierr)
   if (ierr .ne. MPI_SUCCESS) then
      call ESMF_LogWrite( &
           "Problem calling mpi_bcast!", &
           ESMF_LOGMSG_ERROR)
      call endrun(1)
   end if

   use_blacklist = .false.
   if (myid .eq. 0) then
      call esmf_configgetattribute(cf, use_blacklist, &
           label="use_blacklist:", &
           rc=rc)
      if (rc .ne. ESMF_SUCCESS) then
         call ESMF_LogWrite( &
              "Must specify use_blacklist in config file!", &
              ESMF_LOGMSG_ERROR)
         call endrun(1)
      end if
   end if
   call mpi_barrier(mpi_comm_world, ierr)
   call mpi_bcast(use_blacklist, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, &
        ierr)
   if (ierr .ne. MPI_SUCCESS) then
      call ESMF_LogWrite( &
           "Problem calling mpi_bcast!", &
           ESMF_LOGMSG_ERROR)
      call endrun(1)
   end if

   blacklist_file = "NULL"
   if (use_blacklist) then
      if (myid .eq. 0) then
         call esmf_configgetattribute(cf, blacklist_file, &
              label="blacklist_file:", &
              rc=rc)
         if (rc .ne. ESMF_SUCCESS) then
            call ESMF_LogWrite( &
                 "Must specify blacklist_file in config file!", &
                 ESMF_LOGMSG_ERROR)
            call endrun(1)
         end if
      end if
      call mpi_barrier(mpi_comm_world, ierr)
      call mpi_bcast(blacklist_file, len(blacklist_file), MPI_CHARACTER, 0, &
           MPI_COMM_WORLD, ierr)
      if (ierr .ne. MPI_SUCCESS) then
         call ESMF_LogWrite( &
              "Problem calling mpi_bcast!", &
              ESMF_LOGMSG_ERROR)
         call endrun(1)
      end if
   end if

   if (use_blacklist) then
      call fetch_blacklist(blacklist_file, blacklist_stns, nstns)
   end if

   sattype = "NULL"
   if (myid .eq. 0) then
      call esmf_configgetattribute(cf, sattype, &
           label="sattype:", &
           rc=rc)
      if (rc .ne. ESMF_SUCCESS) then
         call ESMF_LogWrite( &
              "Must specify sattype in config file!", &
              ESMF_LOGMSG_ERROR)
         call endrun(1)
      end if
      found_sattype = .false.
      do i = 1, max_sattypes
         if (trim(sattype) .eq. trim(sattypes(i))) then
            found_sattype = .true.
            exit
         end if
      end do
      if (.not. found_sattype) then
         call ESMF_LogWrite( &
              "sattype not supported!", &
              ESMF_LOGMSG_ERROR)
         call ESMF_LogWrite( &
              "sattype must be one of:", &
              ESMF_LOGMSG_ERROR)
         do i = 1, max_sattypes
            call ESMF_LogWrite( trim(sattypes(i)), &
                 ESMF_LOGMSG_ERROR)
         end do
         call endrun(1)
      end if
   end if
   call mpi_barrier(mpi_comm_world, ierr)
   call mpi_bcast(sattype, len(sattype), MPI_CHARACTER, 0, MPI_COMM_WORLD, &
        ierr)
   if (ierr .ne. MPI_SUCCESS) then
      call ESMF_LogWrite( &
           "Problem calling mpi_bcast!", &
           ESMF_LOGMSG_ERROR)
      call endrun(1)
   end if

   imax_lon = 0
   if (myid .eq. 0) then
      call esmf_configgetattribute(cf, imax_lon, &
           label="imax_lon:", &
           rc=rc)
      if (rc .ne. ESMF_SUCCESS) then
         call ESMF_LogWrite( &
              "Must specify imax_lon in config file!", &
              ESMF_LOGMSG_ERROR)
         call endrun(1)
      end if
      if (imax_lon .lt. 1) then
         call ESMF_LogWrite( &
              "Must specify positive imax_lon in config file!", &
              ESMF_LOGMSG_ERROR)
         call endrun(1)
      end if
   end if
   call mpi_barrier(mpi_comm_world, ierr)
   call mpi_bcast(imax_lon, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
   if (ierr .ne. MPI_SUCCESS) then
      call ESMF_LogWrite( &
           "Problem calling mpi_bcast!", &
           ESMF_LOGMSG_ERROR)
      call endrun(1)
   end if

   jmax_lat = 0
   if (myid .eq. 0) then
      call esmf_configgetattribute(cf, jmax_lat, &
           label="jmax_lat:", &
           rc=rc)
      if (rc .ne. ESMF_SUCCESS) then
         call ESMF_LogWrite( &
              "Must specify jmax_lat in config file!", &
              ESMF_LOGMSG_ERROR)
         call endrun(1)
      end if
      if (jmax_lat .lt. 1) then
         call ESMF_LogWrite( &
              "Must specify positive jmax_lat in config file!", &
              ESMF_LOGMSG_ERROR)
         call endrun(1)
      end if
   end if
   call mpi_barrier(mpi_comm_world, ierr)
   call mpi_bcast(jmax_lat, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
   if (ierr .ne. MPI_SUCCESS) then
      call ESMF_LogWrite( &
           "Problem calling mpi_bcast!", &
           ESMF_LOGMSG_ERROR)
      call endrun(1)
   end if

   dlon = 0
   if (myid .eq. 0) then
      call esmf_configgetattribute(cf, dlon, &
           label="dlon:", &
           rc=rc)
      if (rc .ne. ESMF_SUCCESS) then
         call ESMF_LogWrite( &
              "Must specify dlon in config file!", &
              ESMF_LOGMSG_ERROR)
         call endrun(1)
      end if
      if (dlon .lt. 0) then
         call ESMF_LogWrite( &
              "Must specify positive dlon in config file!", &
              ESMF_LOGMSG_ERROR)
         call endrun(1)
      end if
   end if
   call mpi_barrier(mpi_comm_world, ierr)
   call mpi_bcast(dlon, 1, MPI_REAL, 0, MPI_COMM_WORLD, ierr)
   if (ierr .ne. MPI_SUCCESS) then
      call ESMF_LogWrite( &
           "Problem calling mpi_bcast!", &
           ESMF_LOGMSG_ERROR)
      call endrun(1)
   end if

   dlat = 0
   if (myid .eq. 0) then
      call esmf_configgetattribute(cf, dlat, &
           label="dlat:", &
           rc=rc)
      if (rc .ne. ESMF_SUCCESS) then
         call ESMF_LogWrite( &
              "Must specify dlat in config file!", &
              ESMF_LOGMSG_ERROR)
         call endrun(1)
      end if
      if (dlat .lt. 0) then
         call ESMF_LogWrite( &
              "Must specify positive dlat in config file!", &
              ESMF_LOGMSG_ERROR)
         call endrun(1)
      end if
   end if
   call mpi_barrier(mpi_comm_world, ierr)
   call mpi_bcast(dlat, 1, MPI_REAL, 0, MPI_COMM_WORLD, ierr)
   if (ierr .ne. MPI_SUCCESS) then
      call ESMF_LogWrite( &
           "Problem calling mpi_bcast!", &
           ESMF_LOGMSG_ERROR)
      call endrun(1)
   end if

   dist_thresh = 0
   if (myid .eq. 0) then
      call esmf_configgetattribute(cf, dist_thresh, &
           label="dist_thresh:", &
           rc=rc)
      if (rc .ne. ESMF_SUCCESS) then
         call ESMF_LogWrite( &
              "Must specify dist_thresh in config file!", &
              ESMF_LOGMSG_ERROR)
         call endrun(1)
      end if
      if (dist_thresh .lt. 0) then
         call ESMF_LogWrite( &
              "Must specify positive dist_thresh in config file!", &
              ESMF_LOGMSG_ERROR)
         call endrun(1)
      end if
   end if
   call mpi_barrier(mpi_comm_world, ierr)
   call mpi_bcast(dist_thresh, 1, MPI_REAL, 0, MPI_COMM_WORLD, ierr)
   if (ierr .ne. MPI_SUCCESS) then
      call ESMF_LogWrite( &
           "Problem calling mpi_bcast!", &
           ESMF_LOGMSG_ERROR)
      call endrun(1)
   end if

   ! Create ESMF Time and TimeInterval objects
   call esmf_timeset(starttime, yy=startyear, mm=startmonth, dd=startday, &
        h=starthour, m=0, s=0, rc=rc)
   if (rc .ne. ESMF_SUCCESS) then
      call ESMF_LogWrite( &
           "Cannot set starttime object!", &
           ESMF_LOGMSG_ERROR)
      call endrun(1)
   end if
   call esmf_timeset(endtime, yy=endyear, mm=endmonth, dd=endday, &
        h=endhour, m=0, s=0, rc=rc)
   if (rc .ne. ESMF_SUCCESS) then
      call ESMF_LogWrite( &
           "Cannot set endtime object!", &
           ESMF_LOGMSG_ERROR)
      call endrun(1)
   end if

   call esmf_timeintervalset(deltatime, yy=intervalyear, &
        mm=intervalmonth, d=intervalday, h=intervalhour, &
        m=0, s=0, rc=rc)
   if (rc .ne. ESMF_SUCCESS) then
      call ESMF_LogWrite( &
           "Cannot set deltatime object!", &
           ESMF_LOGMSG_ERROR)
      call endrun(1)
   end if

   ! Count the number of possible files
   curtime = starttime
   icount = 0
   do
      if (curtime > endtime) exit
      icount = icount + 1
      curtime = curtime + deltatime
   end do
   max_num_files = icount

   ! Initialize some variables
   iunit_input = 10
   iunit_out_vario = 12

   allocate(icounts_vario(max_vario_bins))
   allocate(vario(max_vario_bins))
   vario(:) = 0
   icounts_vario(:) = 0

   if (adjustl(sattype) == adjustl(SSMI)) then
      call readfiles(myid, numprocs, starttime, endtime, deltatime, datadir, &
           max_stations, max_sat_reports, &
           max_vario_bins, vario_bin_dist, vario, icounts_vario, &
           imax_lon, jmax_lat, dlon, dlat, dist_thresh,&
           is_ssmi, &
           use_blacklist, nstns, blacklist_stns)
   elseif (adjustl(sattype) == adjustl(CMORPH)) then
      call readfiles(myid, numprocs, starttime, endtime, deltatime, datadir, &
           max_stations, max_sat_reports, &
           max_vario_bins, vario_bin_dist, vario, icounts_vario, &
           imax_lon, jmax_lat, dlon, dlat, dist_thresh, &
           is_cmorph, &
           use_blacklist, nstns, blacklist_stns)
   elseif (adjustl(sattype) == adjustl(GEOPRECIP)) then
      call readfiles(myid, numprocs, starttime, endtime, deltatime, datadir, &
           max_stations, max_sat_reports, &
           max_vario_bins, vario_bin_dist, vario, icounts_vario, &
           imax_lon, jmax_lat, dlon, dlat, dist_thresh, &
           is_geoprecip, &
           use_blacklist, nstns, blacklist_stns)
   elseif (adjustl(sattype) == adjustl(IMERG)) then
      call readfiles(myid, numprocs, starttime, endtime, deltatime, datadir, &
           max_stations, max_sat_reports, &
           max_vario_bins, vario_bin_dist, vario, icounts_vario, &
           imax_lon, jmax_lat, dlon, dlat, dist_thresh, &
           is_imerg, &
           use_blacklist, nstns, blacklist_stns)
   else
      call ESMF_LogWrite( &
           'Internal error, cannot search for '//adjustl(sattype), &
           ESMF_LOGMSG_ERROR)
      call endrun(1)
   end if

   if (myid .eq. 0) then
      call ESMF_LogWrite( &
           "Finished reading files", &
           ESMF_LOGMSG_INFO)
   end if

   ! Finish the semivariogram calculation
   if (myid .eq. 0) then
      do j = 1, MAX_VARIO_BINS
         if (icounts_vario(j) .gt. 0) then
            vario(j) = dble(0.5) * vario(j) / dble(icounts_vario(j))
         end if
      end do ! j
   end if

   ! Now write the significant results to file.
   if (myid .eq. 0) then
      open(unit = iunit_out_vario, &
           file = trim(outfile), &
           status="UNKNOWN", &
           iostat = istat)
      if (istat .ne. 0) then
         call ESMF_LogWrite( &
              'Problem opening output file '//trim(outfile), &
              ESMF_LOGMSG_WARNING)
         close(unit=iunit_out_vario)
      end if
      call ESMF_LogWrite( &
           'Writing to output file '//trim(outfile), &
           ESMF_LOGMSG_INFO)
      do j = 1, max_vario_bins
         if (icounts_vario(j) .eq. 0) cycle
         write(iunit_out_vario, '(A,f10.0,A,f8.3,A,I14.14)') &
              ' dist: ', real((j-1)*vario_bin_dist), &
              ' semivariogram: ', vario(j),' icount: ',icounts_vario(j)
      end do ! j
      close(unit=iunit_out_vario)
   end if

   ! The end
   call mpi_barrier(mpi_comm_world, ierr)
   call endrun(0)

contains

   subroutine readfiles(myid,numprocs, starttime, endtime, deltatime, &
        datadir, max_stations, max_sat_reports, max_vario_bins, &
        vario_bin_dist, vario, icounts_vario, &
        imax_lon, jmax_lat, dlon, dlat, dist_thresh, &
        is_sattype, &
        use_blacklist, nstns, blacklist_stns)

      ! Imports
      use USAF_GridHashMod, only: GridHash, newGridHash, destroyGridHash, &
           insertIntoGridHash, getObindexVectorFromGridHash, &
           createIJForGridHash
      use USAF_ReportsMod, only: Reports, newReports, getNobs, getReport, &
           destroyReports, appendToReports, bcast_reports
      use USAF_StationsMod, only: great_circle_distance

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

      ! Local variables
      type(esmf_time) :: curtime
      integer :: yyyy, mm, dd, hh
      character(len=10) :: yyyymmddhh, yyyymmddhh_k, yyyymmddhh_kk
      integer :: rc
      character(len=44) :: infile
      character(len=123) :: fullpath
      integer :: istat
      integer :: iunit_input
      type(Reports) :: R_gages_sat, R_matches
      character(len=80) :: line
      character(len=10) :: platform, platform_k, platform_kk
      character(len=10) :: network, network_k, network_kk, network_new
      real :: latitude, longitude, lat_k, lon_k, lat_kk, lon_kk
      real :: O, B, A, O_k, B_k, O_kk, B_kk, B_new
      real :: min_dist_k_kk, dist_k_kk
      integer :: nobs_gages, nobs_sat
      integer :: nobs
      integer :: i,j,k,kk
      real :: OMB_k, OMB_kk
      integer :: index
      integer*8 :: sampleSize
      integer :: id, id_incr
      double precision, allocatable :: vario_proc(:)
      integer*4, allocatable :: icounts_vario_proc_i4(:)
      double precision, allocatable :: vario_allproc(:)
      integer*4, allocatable :: icounts_vario_allproc_i4(:)
      integer :: sample_size_proc_i4, sample_size_allproc_i4
      integer :: ii,jj
      type(GridHash) :: sat2d, gages2d
      integer, allocatable :: gage_obindexVector(:)
      integer, allocatable :: sat_obindexVector(:)
      integer :: jdelta_lat, idelta_lon
      double precision :: t0, t1, t2, t3
      integer :: icount_gages, icount_sat
      integer :: istn
      integer :: count_skips
      logical :: skip

      iunit_input = 10
      sampleSize = 0

      ! We will loop through each file.  After reading the file, we will
      ! create semivariogram contributions.
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
1000        format(A,I4.4,I2.2,I2.2,I2.2,A)
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
         end if

         ! First, get all the gages and specified satellite obs in the
         ! current file.
         count_skips = 0
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
                           !print*, "[INFO] Skipping blacklisted station ", &
                           !     adjustl(platform)
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
            !print*, "[INFO] Skipped ", count_skips, " blacklisted stations"
            !print*, '[INFO] Working with ', icount_gages, ' gages'
            !print*, '[INFO] Working with ', icount_sat, ' ', trim(sattype), &
            !     ' obs'
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
         sat2d = newGridHash(IMAX_LON, JMAX_LAT, DLON, DLAT)
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

         ! Now we try to match each gage to the nearest sat ob
         call mpi_barrier(mpi_comm_world, ierr)
         R_matches = newReports(MAX_STATIONS)
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
                  end if
               end do ! k
               if (allocated(gage_obindexVector)) then
                  deallocate(gage_obindexVector)
               end if
            end do ! ilon
         end do ! jlat

         ! Clean up
         call destroyReports(R_gages_sat)
         call destroyGridHash(gages2d)
         call destroyGridHash(sat2d)

         ! Update semivariogram sums and counts
         nobs = getNobs(R_matches)
         if (myid .eq. 0) then
            !print*, '[INFO] Working with ', nobs, ' gage/', trim(sattype), &
            !     ' matches...'
         end if
         t0 = mpi_wtime()
         do k = 1, nobs

            allocate(vario_proc(max_vario_bins))
            vario_proc(:) = 0
            allocate(icounts_vario_proc_i4(max_vario_bins))
            icounts_vario_proc_i4(:) = 0
            id = -1
            id_incr = 1
            sample_size_proc_i4 = 0

            call getReport(R_matches, k, platform=platform_k,&
                 latitude=lat_k, longitude=lon_k, &
                 O=O_k, B=B_k)

            OMB_k = O_k - B_k

            t1 = mpi_wtime()

            do kk = k+1, nobs

               ! See if current processor is responsible for this j ob.
               call pick_proc(id, id_incr, numprocs)
               if (id .ne. myid) cycle

               call getReport(R_matches, kk, platform=platform_kk,&
                    latitude=lat_kk, longitude=lon_kk, &
                    O=O_kk, B=B_kk)

               OMB_kk = O_kk - B_kk

               !! Make sure at least one gage has rain.
               !if (.not. (O_i > -3) .or. .not. (O_j > -3)) cycle

               dist_k_kk = great_circle_distance(lat_k, lon_k, lat_kk, lon_kk)
               index = int(dist_k_kk*0.001/vario_bin_dist) + 1
               if (index .gt. max_vario_bins) cycle

               vario_proc(index) = vario_proc(index) + &
                    ((OMB_k-OMB_kk)*(OMB_k-OMB_kk))
               icounts_vario_proc_i4(index) = icounts_vario_proc_i4(index) + 1
               sample_size_proc_i4 = sample_size_proc_i4 + 1
            end do ! i

            ! Now need to collect the data from each processor
            sample_size_allproc_i4 = 0
            call mpi_barrier(mpi_comm_world, ierr)
            call mpi_allreduce(sample_size_proc_i4, sample_size_allproc_i4, &
                 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
            sampleSize = sampleSize + sample_size_allproc_i4

            allocate(vario_allproc(max_vario_bins))
            vario_allproc(:) = 0
            call mpi_barrier(mpi_comm_world,ierr)
            call mpi_allreduce(vario_proc, vario_allproc, max_vario_bins, &
                 MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
            vario(:) = vario(:) + vario_allproc(:)
            deallocate(vario_proc)
            deallocate(vario_allproc)

            allocate(icounts_vario_allproc_i4(max_vario_bins))
            icounts_vario_allproc_i4(:) = 0
            call mpi_barrier(mpi_comm_world, ierr)
            call mpi_allreduce(icounts_vario_proc_i4, &
                 icounts_vario_allproc_i4, &
                 max_vario_bins, &
                 MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
            icounts_vario(:) = icounts_vario(:) + icounts_vario_allproc_i4(:)
            deallocate(icounts_vario_proc_i4)
            deallocate(icounts_vario_allproc_i4)

            t2 = mpi_wtime()
            !print*,'Ob j  ',j,' took ',t2-t1,' seconds'

         end do ! k

         t3 = mpi_wtime()
         !print*, '[INFO] Processing completed after ', t3-t0, ' seconds'

         ! Clean up for next file
         call destroyReports(R_matches)

         if (myid .eq. 0) then
            !print*, '[INFO] Semivariogram sample size now ', sampleSize
         end if

       enddo ! next file
    end subroutine readfiles

end program main
