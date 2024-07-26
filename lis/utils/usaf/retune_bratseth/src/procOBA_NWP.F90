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
! ROUTINE: procOBA_NWP
!
! DESCRIPTION:
! Main program for generating semivariogram of observation-NWP differences
! with distance.
!
! REVISION HISTORY:
! 26 0ct 2020:  Eric Kemp.  Initial Specification. Multi-process
!   MPI doesn't work yet, just use single process.
! 15 Dec 2020:  Eric Kemp.  Added user-defined logfile name.

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
   integer, parameter :: maxlen_obtype = 6
   character(len=maxlen_obtype) :: obtype
   integer, parameter :: max_obtypes = 4
   character(len=maxlen_obtype) :: obtypes(max_obtypes)
   logical :: found_obtype
   double precision :: t0, t1, t2, t3
   integer :: numprocs, ierr
   integer :: myid
   logical :: is_precip
   integer :: num_args
   logical :: use_blacklist
   character(len=255) :: blacklist_file
   character(len=9), allocatable :: blacklist_stns(:)
   integer :: nstns
   character(len=255) :: logname

   ! TODO:  Add other ob types.
   character(len=maxlen_obtype), parameter :: T2M     = "T2M"
   character(len=maxlen_obtype), parameter :: RH2M    = "RH2M"
   character(len=maxlen_obtype), parameter :: SPD10M  = "SPD10M"
   character(len=maxlen_obtype), parameter :: GAGE   =  "Gage"
   obtypes = (/T2M, RH2M, SPD10M, GAGE/)

   ! Set logfile name
   logname = 'procOBA_NWP.log'
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
      call esmf_ConfigGetAttribute(cf, max_stations, label="max_stations:", &
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

   max_vario_bins = 0
   if (myid .eq. 0) then
      call esmf_ConfigGetAttribute(cf, max_vario_bins, &
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

   obtype = "NULL"
   if (myid .eq. 0) then
      call esmf_configgetattribute(cf, obtype, &
           label="obtype:", &
           rc=rc)
      if (rc .ne. ESMF_SUCCESS) then
         call ESMF_LogWrite( &
              "Must specify obtype in config file!", &
              ESMF_LOGMSG_ERROR)
         call endrun(1)
      end if
      found_obtype = .false.
      do i = 1, max_obtypes
         if (trim(obtype) .eq. trim(obtypes(i))) then
            found_obtype = .true.
            exit
         end if
      end do
      if (.not. found_obtype) then
         call ESMF_LogWrite( &
              "obtype not supported!", &
              ESMF_LOGMSG_ERROR)
         call ESMF_LogWrite( &
              "obtype must be one of:", &
              ESMF_LOGMSG_ERROR)
         do i = 1, max_obtypes
            call ESMF_LogWrite( &
                 trim(obtypes(i)), &
                 ESMF_LOGMSG_ERROR)
         end do
         call endrun(1)
      end if
   end if
   call mpi_barrier(mpi_comm_world, ierr)
   call mpi_bcast(obtype, len(obtype), MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
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

   ! TODO:  Add support for other observations
   if ( adjustl(obtype) == adjustl(T2M) .or. &
        adjustl(obtype) == adjustl(RH2M) .or. &
        adjustl(obtype) == adjustl(SPD10M)) then
      is_precip = .false.
      call readfiles(myid, numprocs, starttime, endtime, deltatime, datadir, &
           max_stations, &
           max_vario_bins, vario_bin_dist, vario, icounts_vario, &
           is_precip, is_gage, &
           use_blacklist, nstns, blacklist_stns)
   else if (adjustl(obtype) == adjustl(GAGE)) then
      is_precip = .true.
      call readfiles(myid, numprocs, starttime, endtime, deltatime, datadir, &
           max_stations, &
           max_vario_bins, vario_bin_dist, vario, icounts_vario, &
           is_precip, is_gage, &
           use_blacklist, nstns, blacklist_stns)
   else
      call ESMF_LogWrite( &
           "Internal error, cannot search for "//adjustl(obtype), &
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
      open(unit=iunit_out_vario, &
           file=trim(outfile), &
           status="UNKNOWN", &
           iostat=istat)
      if (istat .ne. 0) then
         call ESMF_LogWrite( &
              "Problem opening output file "//trim(outfile), &
              ESMF_LOGMSG_WARNING)
         close(unit=iunit_out_vario)
      end if
      call ESMF_LogWrite( &
           "Writing to output file "//trim(outfile), &
           ESMF_LOGMSG_INFO)
      do j = 1, max_vario_bins
         if (icounts_vario(j) .eq. 0) cycle
         write(iunit_out_vario, '(A,f10.0,A,f7.3,A,I14.14)') &
              ' dist: ',  real((j-1)*vario_bin_dist), &
              ' semivariogram: ', vario(j), ' icount: ', icounts_vario(j)
      end do ! j
      close(unit=iunit_out_vario)
   end if

   ! The end
   call mpi_barrier(mpi_comm_world, ierr)
   call endrun(0)

contains

   subroutine readfiles(myid, numprocs, starttime, endtime, deltatime, &
        datadir, max_stations, max_vario_bins, vario_bin_dist, vario, &
        icounts_vario, is_precip, is_obtype, &
        use_blacklist, nstns, blacklist_stns)

      ! Imports
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
      integer, intent(in) :: max_vario_bins
      integer, intent(in) :: vario_bin_dist
      double precision, intent(inout) :: vario(max_vario_bins)
      integer*8, intent(inout) :: icounts_vario(max_vario_bins)
      logical, intent(in) :: is_precip
      logical, external :: is_obtype ! Function
      logical, intent(in) :: use_blacklist
      integer, intent(in) :: nstns
      character(len=9), allocatable, intent(in) ::  blacklist_stns(:)

      ! Local variables
      type(esmf_time) :: curtime
      integer :: yyyy, mm, dd, hh
      character(len=10) :: yyyymmddhh
      integer :: rc
      character(len=44) :: infile
      character(len=123) :: fullpath
      integer :: istat
      integer :: iunit_input
      type(Reports) :: R_obs
      character(len=80) :: line
      character(len=10) :: platform, platform_i,platform_j
      character(len=10) :: network
      real :: latitude, longitude, lat_i, lon_i, lat_j, lon_j
      real :: O, B, A, O_i, B_i, O_j, B_j
      integer :: nobs
      integer :: i,j
      real :: OMB_i, OMB_j
      real :: dist_ij
      integer :: index
      integer*8 :: sampleSize
      integer :: id, id_incr
      double precision, allocatable :: vario_proc(:)
      integer*4, allocatable :: icounts_vario_proc_i4(:)
      double precision, allocatable :: vario_allproc(:)
      integer*4, allocatable :: icounts_vario_allproc_i4(:)
      integer :: sample_size_proc_i4, sample_size_allproc_i4
      integer :: istn
      integer :: count_skips
      logical :: skip

      iunit_input = 10
      sampleSize = 0

      ! We will loop through each file.  After reading the file, we will
      ! create semivariogram contributions.  By doing this for each file, we
      ! will hopefully reduce the number of obs that must be compared for a
      ! particular date and time.
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
            if (is_precip) then
               write(infile,1000) 'oba_', yyyy, mm, dd, hh, '_12.txt'
            else
               write(infile,1000) 'oba_', yyyy, mm, dd, hh, '_01.txt'
            end if
1000        format(A, I4.4, I2.2, I2.2, I2.2, A)
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
                 "Problem opening "//trim(fullpath), &
                 ESMF_LOGMSG_WARNING)
               close(unit=iunit_input)
            end if
            cycle
         end if

         if (myid .eq. 0) then
            call ESMF_LogWrite( &
                 "Reading file "//trim(infile), &
                 ESMF_LOGMSG_INFO)
         end if

         ! First, get all the data in the current file.
         count_skips = 0
         R_obs = newReports(max_stations)
         if (myid .eq. 0) then
            do
               read(iunit_input, '(A)', iostat=istat) line
               if (istat .ne. 0) exit
               if (line(2:2) .eq. '#') cycle
               read(line, &
                    '(A10,1x,A10,1x,F8.3,1x,F8.3,1x,F8.3,1x,F8.3,1x,F8.3)') &
                    network, platform, latitude, longitude, O, B, A

               ! Only save requested obtype
               if (.not. is_obtype(adjustl(network))) cycle

               ! Apply blacklist if requested
               skip = .false.
               if (use_blacklist) then
                  do istn = 1, nstns
                     !print*, "EMK: blacklist, platform = ", &
                     !     trim(blacklist_stns(istn)), " ", adjustl(platform)
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

               call appendToReports(R_obs, network, platform, yyyymmddhh, &
                    latitude, longitude, O, B, A)
            end do ! loop through file
            close(unit=iunit_input)
         end if

         !print*, "[INFO] Skipped ", count_skips, " blacklisted stations"

         ! Share reports across processors
         call mpi_barrier(mpi_comm_world, ierr)
         call bcast_reports(R_obs, myid, ierr)
         if (ierr .ne. MPI_SUCCESS) then
            call ESMF_LogWrite( &
                 "Cannot share reports across processors!", &
                 ESMF_LOGMSG_ERROR)
            call endrun(1)
         end if

         ! Update semivariogram sums and counts
         nobs = getNobs(R_obs)
         if (myid .eq. 0) then
            !print*, '[INFO] Working with ', nobs, ' ', trim(obtype), &
            !     '/NWP matches...'
         end if
         t0 = mpi_wtime()
         do j = 1,nobs

            allocate(vario_proc(max_vario_bins))
            vario_proc(:) = 0
            allocate(icounts_vario_proc_i4(max_vario_bins))
            icounts_vario_proc_i4(:) = 0
            id = -1
            id_incr = 1
            sample_size_proc_i4 = 0

            call getReport(R_obs, j, platform=platform_j,&
                 latitude=lat_j, longitude=lon_j, &
                 O=O_j, B=B_j)

            OMB_j = O_j - B_j

            t1 = mpi_wtime()

            do i = j+1,nobs

               ! See if current processor is responsible for this j ob.
               call pick_proc(id,id_incr, numprocs)
               if (id .ne. myid) cycle

               call getReport(R_obs,i, platform=platform_i,&
                    latitude=lat_i, longitude=lon_i, &
                    O=O_i, B=B_i)

               OMB_i = O_i - B_i

               !! Make sure at least one gage has rain.
               !if (.not. (O_i > -3) .or. .not. (O_j > -3)) cycle

               dist_ij = great_circle_distance(lat_i, lon_i, lat_j, lon_j)
               index = int(dist_ij*0.001/vario_bin_dist) + 1
               if (index .gt. max_vario_bins) cycle

               vario_proc(index) = vario_proc(index) + &
                    ((OMB_i-OMB_j)*(OMB_i-OMB_j))
               icounts_vario_proc_i4(index) = icounts_vario_proc_i4(index) + 1
               sample_size_proc_i4 = sample_size_proc_i4 + 1
            end do ! i

            ! Now need to collect the data from each processor
            sample_size_allproc_i4 = 0
            call mpi_barrier(mpi_comm_world,ierr)
            call mpi_allreduce(sample_size_proc_i4, sample_size_allproc_i4, &
                 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
            sampleSize = sampleSize + sample_size_allproc_i4

            allocate(vario_allproc(max_vario_bins))
            vario_allproc(:) = 0
            call mpi_barrier(mpi_comm_world,ierr)
            call mpi_allreduce(vario_proc,vario_allproc, max_vario_bins, &
                 MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
            vario(:) = vario(:) + vario_allproc(:)
            deallocate(vario_proc)
            deallocate(vario_allproc)

            allocate(icounts_vario_allproc_i4(max_vario_bins))
            icounts_vario_allproc_i4(:) = 0
            call mpi_barrier(mpi_comm_world,ierr)
            call mpi_allreduce(icounts_vario_proc_i4, &
                 icounts_vario_allproc_i4, &
                 max_vario_bins, &
                 MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
            icounts_vario(:) = icounts_vario(:) + icounts_vario_allproc_i4(:)
            deallocate(icounts_vario_proc_i4)
            deallocate(icounts_vario_allproc_i4)

            t2 = mpi_wtime()
            !print*,'Ob j  ',j,' took ',t2-t1,' seconds'

         end do ! j

         t3 = mpi_wtime()
         !print*,'[INFO] Processing completed after ', t3-t0, ' seconds'

         ! Clean up for next file
         call destroyReports(R_obs)

         if (myid .eq. 0) then
            !print*, '[INFO] Semivariogram sample size now ', sampleSize
         end if
      end do

   end subroutine readfiles

end program main
