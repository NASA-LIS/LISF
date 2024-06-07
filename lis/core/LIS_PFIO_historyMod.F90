!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#if ( defined USE_PFIO )
#include "MAPL_ErrLog.h"
#include "unused_dummy.H"
#endif

#include "LIS_misc.h"
module LIS_PFIO_historyMod
!BOP
!
! !MODULE: LIS_PFIO_historyMod
! 
! !DESCRIPTION: 
!  Provides interfaces to manage LIS model output 
!  (LIS History) using GEOS MAPL/PFIO. 
!  The output files are in netCDF format.
!  
! !USES: 
   use LIS_coreMod
   use LIS_histDataMod
   use LIS_timeMgrMod
   use LIS_logMod

   use LIS_mpiMod
#if ( defined USE_PFIO )
   use, intrinsic :: iso_fortran_env, only: REAL32
   use MAPL
   use ESMF
   use pFIO_ClientManagerMod, only: o_Clients
   use pFIO_UnlimitedEntityMod
   use LIS_PFIO_utilsMod, only : PFIO_write_var

   implicit none

   PRIVATE
!--------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!--------------------------------
   public :: LIS_PFIO_create_file_metadata
   public :: LIS_PFIO_write_data
   public :: LIS_PFIO_create_routing_metadata
   public :: LIS_PFIO_write_routingdata

   REAL               ::          pfio_vmin 
   REAL               ::          pfio_vmax
   REAL               :: pfio_missing_value
   REAL               ::    pfio_fill_value
   REAL, DIMENSION(2) ::   pfio_valid_range 


   INTEGER, PARAMETER :: MAX_NUM_VLEVELS = 12
   integer            :: idx_field2d        ! index of 2d field to ne written out
   integer            :: idx_field3d        ! index of 3d field to ne written out
   integer            :: idx_field4d        ! index of 4d field to ne written out
   INTEGER            :: total_num_2Dfields(10) ! total number of 2D fields to be written out
   INTEGER            :: total_num_3Dfields(10) ! total number of 3D fields to be written out
   INTEGER            :: total_num_4Dfields(10) ! total number of 4D fields to be written out

   ! Derived type to hold the values of the fields to be written

   TYPE local_2Dfield_var
      real, allocatable :: var2d(:,:)
      integer           :: idx_field2d
   END TYPE local_2Dfield_var

   TYPE local_3Dfield_var
      real, allocatable :: var3d(:,:,:)
      integer           :: idx_field3d
   END TYPE local_3Dfield_var

   TYPE local_4Dfield_var
      real, allocatable :: var4d(:,:,:,:)
      integer           :: idx_field4d
   END TYPE local_4Dfield_var

   type(LIS_metadataEntry), pointer :: xlat, xlong

#else
   integer :: dummy_int
#endif

#if ( defined USE_PFIO )
!EOP
!------------------------------------------------------------------------------ 
CONTAINS
!------------------------------------------------------------------------------ 
!BOP
! !ROUTINE: LIS_PFIO_create_file_metadata
! \label{LIS_PFIO_create_file_metadata}
! 
! !INTERFACE: LIS_PFIO_create_file_metadata
   subroutine LIS_PFIO_create_file_metadata(n, PFIOmodel_idx, outInterval, &
         nsoillayers, lyrthk, model_name, group)
! !USES: 
!
! !INPUTS PARAMETERS: 
      integer,          intent(in) :: n 
      integer,          intent(in) :: PFIOmodel_idx 
      real,             intent(in) :: outInterval
      integer,          intent(in) :: nsoillayers
      real,             intent(in) :: lyrthk(nsoillayers)
      character(len=*), intent(in), optional :: model_name
      integer,          intent(in), optional :: group
! 
! !DESCRIPTION: 
!  This routine creates a PFIO object that contains metadata to be used
!  in a netCDF file.
!  The arguments are: 
!  \begin{description}
!    \item[n] index of the nest
!    \item[ftn] file unit for the output file
!    \item[ftn\_stats] file unit for the output statistics file
!    \item[outInterval]   history output frequency
!    \item[nsoillayers]  Number of soil layers
!    \item[lyrthk]   Thickness of soil layers
!    \item[model\_name] Name of the model that generates the output
!  \end{description}
!
!   The routines invoked are: 
!   \begin{description}
!   \item[defineNETCDFheadervar](\ref{PFIO_define_variable_header})
!     writes the required headers for a single variable
!   \item[LIS\_verify](\ref{LIS_verify})
!     call to check if the return value is valid or not.
!   \end{description}
!
! !LOCAL VARIABLES:
      Type(Variable)          :: v

      integer                 :: dimID(4)
      integer                 :: tdimID,xtimeID,ensID
      integer                 :: t,c,r,i,index1
      real, allocatable       :: ensval(:) 
      character(len=8)        :: xtime_begin_date
      character(len=6)        :: xtime_begin_time
      character(len=50)       :: xtime_units
      character(len=50)       :: xtime_timeInc
      integer                 :: iret
      integer                 :: group_
      character(len=100)      :: model_name_
      integer                 :: status, gindex
      ! Note that the fix to add lat/lon to the NETCDF output will output
      ! undefined values for the water points. 
      character(len=8)        :: date
      character(len=10)       :: time
      character(len=90)       :: name_dims(4)
      character(len=5)        :: zone
      integer, dimension(8)   :: values
      type(LIS_metadataEntry), pointer :: dataEntry
      real                    :: time_data(1)
      REAL                    :: SOUTH_WEST_CORNER_LAT
      REAL                    :: SOUTH_WEST_CORNER_LON
      integer                 :: ierr
      integer                 :: vcol_id
!EOP
!---------------------------------------------------------------------------------------------
!BOC
      write(LIS_logunit,'(a,i2,a)')'[INFO-PFIO] Creating the file metadata [',PFIOmodel_idx,']'
      !--> ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      !--> We need to STOP if the "1d tilespace" configuration is selected.
      !--> PFIO cannot handle the way local tiles are mapped into the global one.
      !--> ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      IF (LIS_rc%wopt.eq."1d tilespace") THEN
         ierr = 1
         call LIS_verify(ierr,'The 1d tilespace configuration is not supported with PFIO!')
      ENDIF

      if (.NOT.PRESENT(group)) then
         group_ = 1
      else
         group_ = group
      endif

      if (.NOT.PRESENT(model_name)) then
         model_name_ = "model_not_specified"
      else
         model_name_ = model_name
      endif

      call LIS_rescaleCount(n, group_)

      call date_and_time(date,time,zone,values)

      pfio_vmin          = 0.0 ! -MAPL_UNDEF
      pfio_vmax          = 0.0 !  MAPL_UNDEF
      pfio_missing_value = LIS_rc%udef ! MAPL_UNDEF
      pfio_fill_value    = LIS_rc%udef ! MAPL_UNDEF
      pfio_valid_range   = (/-MAPL_UNDEF, MAPL_UNDEF/)

      !----------------------------
      ! Longitude and Latitude data
      !----------------------------

      allocate(xlat)
      allocate(xlong)

      xlat%short_name = "lat"
      xlat%long_name = "latitude"
      xlat%standard_name = "latitude"
      xlat%units = "degree_north"
      xlat%nunits = 1
      xlat%format = 'F'
      xlat%form = 1
      xlat%vlevels = 1
      xlat%timeAvgOpt = 0
      xlat%selectOpt = 1
      xlat%minMaxOpt = 0
      xlat%stdOpt = 0
      allocate(xlat%modelOutput(1,LIS_rc%ntiles(n),xlat%vlevels))
      allocate(xlat%count(1,xlat%vlevels))
      xlat%count = 1
      allocate(xlat%unittypes(1))
      xlat%unittypes(1) = "degree_north"
      xlat%valid_min = 0.0
      xlat%valid_max = 0.0

      xlong%short_name = "lon"
      xlong%long_name = "longitude"
      xlong%standard_name = "longitude"
      xlong%units = "degree_east"
      xlong%nunits = 1
      xlong%format = 'F'
      xlong%form = 1
      xlong%vlevels = 1
      xlong%timeAvgOpt = 0
      xlong%selectOpt = 1
      xlong%minMaxOpt = 0
      xlong%stdOpt = 0
      allocate(xlong%modelOutput(1,LIS_rc%ntiles(n),xlong%vlevels))
      allocate(xlong%count(1,xlong%vlevels))
      xlong%count = 1
      allocate(xlong%unittypes(1))
      xlong%unittypes(1) = "degree_east"
      xlong%valid_min = 0.0
      xlong%valid_max = 0.0

      do t=1,LIS_rc%ntiles(n)
         c = LIS_domain(n)%tile(t)%col
         r = LIS_domain(n)%tile(t)%row
         index1 = LIS_domain(n)%gindex(c,r)
         xlat%modelOutput(1,t,1)  = LIS_domain(n)%grid(index1)%lat
         xlong%modelOutput(1,t,1) = LIS_domain(n)%grid(index1)%lon
      enddo

      ! Write variable data (both non-model and model-based):
      if (LIS_rc%wopt.eq."2d ensemble gridspace") then
         allocate(ensval(LIS_rc%nensem(n)))
         do i = 1, LIS_rc%nensem(n)
            ensval(i) = float(i)
         end do
      endif

      ! defining time field
      write(xtime_units, 200) LIS_rc%yr, LIS_rc%mo, LIS_rc%da, LIS_rc%hr, LIS_rc%mn, LIS_rc%ss
      200   format ('minutes since ',I4.4,'-',I2.2,'-',I2.2,' ',I2.2,':', I2.2,':',I2.2)
      write(xtime_begin_date, fmt='(I4.4,I2.2,I2.2)') LIS_rc%yr, LIS_rc%mo, LIS_rc%da
      write(xtime_begin_time, fmt='(I2.2,I2.2,I2.2)') LIS_rc%hr, LIS_rc%mn, LIS_rc%ss
      write(xtime_timeInc, fmt='(I20)')  nint(outInterval)

      COL_LOOPS: DO vcol_id = 1, LIS_rc%n_vcollections
         write(LIS_logunit,'(a,i2,a1,i2,a2,i2,a)')'[INFO-PFIO] Create file metadata of collection ', vcol_id, '/', LIS_rc%n_vcollections,' [',PFIOmodel_idx,']'
         name_dims(:) = ''

         !------------------
         ! Define dimensions
         !------------------
         if (LIS_rc%wopt.eq."1d tilespace") then 
            name_dims(1) = 'ntiles'
            call PFIO_bundle%fmd(n,vcol_id,PFIOmodel_idx)%add_dimension(TRIM(name_dims(1)), LIS_rc%glbntiles(n), rc=status)
         elseif ( (LIS_rc%wopt.eq."2d gridspace") .OR. (LIS_rc%wopt.eq."2d ensemble gridspace") )then 
            name_dims(1) = 'east_west'
            name_dims(2) = 'north_south'
            call PFIO_bundle%fmd(n,vcol_id,PFIOmodel_idx)%add_dimension(TRIM(name_dims(1)), LIS_rc%gnc(n), rc=status) 
            call PFIO_bundle%fmd(n,vcol_id,PFIOmodel_idx)%add_dimension(TRIM(name_dims(2)), LIS_rc%gnr(n), rc=status)

            if (LIS_rc%wopt.eq."2d ensemble gridspace") then 
               name_dims(3) = 'ensemble'
               call PFIO_bundle%fmd(n,vcol_id,PFIOmodel_idx)%add_dimension(TRIM(name_dims(3)), LIS_rc%nensem(n), rc=status)

               v = Variable(type=PFIO_REAL32, dimensions='ensemble')
               call v%add_attribute("units", "ensemble number")
               call v%add_attribute("long_name", "Ensemble numbers")
               call v%add_const_value(UnlimitedEntity(ensval))
               call PFIO_bundle%fmd(n,vcol_id,PFIOmodel_idx)%add_variable('ensemble', v)
            endif
         endif

         ! LIS output is always writing output for a single time record
         call PFIO_bundle%fmd(n,vcol_id,PFIOmodel_idx)%add_dimension('time', pFIO_UNLIMITED, rc=status)

         !------------------------------------------
         ! Define variables and variables attributes
         !------------------------------------------

         ! Time: This is the initial step. 
         !       The attributes will be overwritten each time a file is created.
         time_data = 0.0
         v = Variable(type=PFIO_REAL32, dimensions='time')
         call v%add_attribute("units", trim(xtime_units))
         call v%add_attribute("long_name", "time")
         call v%add_attribute("time_increment", trim(adjustl(xtime_timeInc)))
         call v%add_attribute("begin_date", xtime_begin_date)
         call v%add_attribute("begin_time", xtime_begin_time)
         !call v%add_const_value(UnlimitedEntity( time_data ))
         call PFIO_bundle%fmd(n,vcol_id,PFIOmodel_idx)%add_variable('time', v)

         ! Pointer to header information
         SELECT CASE (group_)
         CASE(1) ! LSM output
            dataEntry => LIS_histData(n)%head_lsm_list
            write(LIS_logunit,'(a,i2,a)')'[INFO-PFIO] Doing HEAD_LSM_LIST [',PFIOmodel_idx,']'
         CASE(2) ! ROUTING
            dataEntry => LIS_histData(n)%head_routing_list
            write(LIS_logunit,'(a,i2,a)')'[INFO-PFIO] Doing HEAD_ROUTING_LIST [',PFIOmodel_idx,']'
         CASE(3) ! RTM
            dataEntry => LIS_histData(n)%head_rtm_list
            write(LIS_logunit,'(a,i2,a)')'[INFO-PFIO] Doing HEAD_RTM_LIST [',PFIOmodel_idx,']'
         CASE(4) ! Irrigation
            dataEntry => LIS_histData(n)%head_irrig_list
            write(LIS_logunit,'(a,i2,a)')'[INFO-PFIO] Doing HEAD_IRRIG_LIST [',PFIOmodel_idx,']'
         END SELECT

         total_num_2Dfields(PFIOmodel_idx) = 0
         total_num_3Dfields(PFIOmodel_idx) = 0
         total_num_4Dfields(PFIOmodel_idx) = 0

         !if ( (LIS_rc%wopt.eq."2d gridspace") .OR. (LIS_rc%wopt.eq."2d ensemble gridspace") )then
         call PFIO_define_variable_header(n, PFIOmodel_idx, vcol_id, xlat,  name_dims, non_model_fields = 1)
         call PFIO_define_variable_header(n, PFIOmodel_idx, vcol_id, xlong, name_dims, non_model_fields = 2)
         !endif

         do while ( associated(dataEntry) )
            call PFIO_define_variable_header(n, PFIOmodel_idx, vcol_id, dataEntry, name_dims)
            dataEntry => dataEntry%next
         enddo

         !------------------
         ! Global attributes
         !------------------
         call PFIO_bundle%fmd(n,vcol_id,PFIOmodel_idx)%add_attribute("missing_value", pfio_missing_value) ! LIS_rc%udef)
         call PFIO_bundle%fmd(n,vcol_id,PFIOmodel_idx)%add_attribute("NUM_SOIL_LAYERS", nsoillayers)
         call PFIO_bundle%fmd(n,vcol_id,PFIOmodel_idx)%add_attribute("SOIL_LAYER_THICKNESSES", lyrthk)
         call PFIO_bundle%fmd(n,vcol_id,PFIOmodel_idx)%add_attribute("title", "LIS land surface model output")
         call PFIO_bundle%fmd(n,vcol_id,PFIOmodel_idx)%add_attribute("institution", trim(LIS_rc%institution))
         call PFIO_bundle%fmd(n,vcol_id,PFIOmodel_idx)%add_attribute("source", trim(model_name_))
         call PFIO_bundle%fmd(n,vcol_id,PFIOmodel_idx)%add_attribute("history",  & 
            "created on date: "//date(1:4)//"-"//date(5:6)//"-"// date(7:8)//"T"//time(1:2)//":"//time(3:4)//":"//time(5:10))
         call PFIO_bundle%fmd(n,vcol_id,PFIOmodel_idx)%add_attribute("references", "Kumar_etal_EMS_2006, Peters-Lidard_etal_ISSE_2007")
         call PFIO_bundle%fmd(n,vcol_id,PFIOmodel_idx)%add_attribute("conventions", "CF-1.6")
         call PFIO_bundle%fmd(n,vcol_id,PFIOmodel_idx)%add_attribute("comment", "website: http://lis.gsfc.nasa.gov/")

         ! Grid information

         !============================================================================
         ! Make sure that all the processes have the same value the lower left corner.
         ! Broadcast the the root process values to all the other processes.
         ! This is only needed to write global attributes in the file.
         !============================================================================
         IF (LIS_masterproc) THEN
            SOUTH_WEST_CORNER_LAT = LIS_rc%gridDesc(n,4)
            SOUTH_WEST_CORNER_LON = LIS_rc%gridDesc(n,5)
         ELSE
            SOUTH_WEST_CORNER_LAT = -9999.0
            SOUTH_WEST_CORNER_LON = -9999.0
         ENDIF
         call MPI_Bcast(SOUTH_WEST_CORNER_LAT,1, MPI_REAL, 0, LIS_mpi_comm, ierr)
         call MPI_Bcast(SOUTH_WEST_CORNER_LON,1, MPI_REAL, 0, LIS_mpi_comm, ierr)

         if (LIS_rc%lis_map_proj.eq."latlon") then   ! latlon
            call PFIO_bundle%fmd(n,vcol_id,PFIOmodel_idx)%add_attribute("MAP_PROJECTION", "EQUIDISTANT CYLINDRICAL")
            call PFIO_bundle%fmd(n,vcol_id,PFIOmodel_idx)%add_attribute("SOUTH_WEST_CORNER_LAT", SOUTH_WEST_CORNER_LAT)
            call PFIO_bundle%fmd(n,vcol_id,PFIOmodel_idx)%add_attribute("SOUTH_WEST_CORNER_LON", SOUTH_WEST_CORNER_LON)
            call PFIO_bundle%fmd(n,vcol_id,PFIOmodel_idx)%add_attribute("DX", LIS_rc%gridDesc(n,9))
            call PFIO_bundle%fmd(n,vcol_id,PFIOmodel_idx)%add_attribute("DY", LIS_rc%gridDesc(n,10))       
         elseif (LIS_rc%lis_map_proj.eq."mercator") then 
            call PFIO_bundle%fmd(n,vcol_id,PFIOmodel_idx)%add_attribute("MAP_PROJECTION", "MERCATOR")
            call PFIO_bundle%fmd(n,vcol_id,PFIOmodel_idx)%add_attribute("SOUTH_WEST_CORNER_LAT", SOUTH_WEST_CORNER_LAT)
            call PFIO_bundle%fmd(n,vcol_id,PFIOmodel_idx)%add_attribute("SOUTH_WEST_CORNER_LON", SOUTH_WEST_CORNER_LON)
            call PFIO_bundle%fmd(n,vcol_id,PFIOmodel_idx)%add_attribute("TRUELAT1", LIS_rc%gridDesc(n,10))
            call PFIO_bundle%fmd(n,vcol_id,PFIOmodel_idx)%add_attribute("STANDARD_LON", LIS_rc%gridDesc(n,11))
            call PFIO_bundle%fmd(n,vcol_id,PFIOmodel_idx)%add_attribute("DX", LIS_rc%gridDesc(n,8))
            call PFIO_bundle%fmd(n,vcol_id,PFIOmodel_idx)%add_attribute("DY", LIS_rc%gridDesc(n,9))
         elseif (LIS_rc%lis_map_proj.eq."lambert") then ! lambert conformal
            call PFIO_bundle%fmd(n,vcol_id,PFIOmodel_idx)%add_attribute("MAP_PROJECTION",  "LAMBERT CONFORMAL")
            call PFIO_bundle%fmd(n,vcol_id,PFIOmodel_idx)%add_attribute("SOUTH_WEST_CORNER_LAT", SOUTH_WEST_CORNER_LAT)
            call PFIO_bundle%fmd(n,vcol_id,PFIOmodel_idx)%add_attribute("SOUTH_WEST_CORNER_LON", SOUTH_WEST_CORNER_LON)
            call PFIO_bundle%fmd(n,vcol_id,PFIOmodel_idx)%add_attribute("TRUELAT1", LIS_rc%gridDesc(n,10))
            call PFIO_bundle%fmd(n,vcol_id,PFIOmodel_idx)%add_attribute("TRUELAT2", LIS_rc%gridDesc(n,7))
            call PFIO_bundle%fmd(n,vcol_id,PFIOmodel_idx)%add_attribute("STANDARD_LON", LIS_rc%gridDesc(n,11))
            call PFIO_bundle%fmd(n,vcol_id,PFIOmodel_idx)%add_attribute("DX", LIS_rc%gridDesc(n,8))
            call PFIO_bundle%fmd(n,vcol_id,PFIOmodel_idx)%add_attribute("DY", LIS_rc%gridDesc(n,9))
         elseif (LIS_rc%lis_map_proj.eq."polar") then ! polar stereographic
            call PFIO_bundle%fmd(n,vcol_id,PFIOmodel_idx)%add_attribute("MAP_PROJECTION", "POLAR STEREOGRAPHIC")
            call PFIO_bundle%fmd(n,vcol_id,PFIOmodel_idx)%add_attribute("SOUTH_WEST_CORNER_LAT", SOUTH_WEST_CORNER_LAT)
            call PFIO_bundle%fmd(n,vcol_id,PFIOmodel_idx)%add_attribute("SOUTH_WEST_CORNER_LON", SOUTH_WEST_CORNER_LON)
            call PFIO_bundle%fmd(n,vcol_id,PFIOmodel_idx)%add_attribute("TRUELAT1", LIS_rc%gridDesc(n,10))
            call PFIO_bundle%fmd(n,vcol_id,PFIOmodel_idx)%add_attribute("ORIENT", LIS_rc%gridDesc(n,7))
            call PFIO_bundle%fmd(n,vcol_id,PFIOmodel_idx)%add_attribute("STANDARD_LON", LIS_rc%gridDesc(n,11))
            call PFIO_bundle%fmd(n,vcol_id,PFIOmodel_idx)%add_attribute("DX", LIS_rc%gridDesc(n,8))
            call PFIO_bundle%fmd(n,vcol_id,PFIOmodel_idx)%add_attribute("DY", LIS_rc%gridDesc(n,9))
         endif       

         ! Create the history collection
         PFIO_bundle%hist_id(n,vcol_id,PFIOmodel_idx) = o_Clients%add_hist_collection(PFIO_bundle%fmd(n,vcol_id,PFIOmodel_idx))
      ENDDO COL_LOOPS

      ! ----> Done with the creation of the PFIO object

      if ( LIS_rc%wopt.eq."2d ensemble gridspace" ) then
         deallocate(ensval)
      endif

   end subroutine LIS_PFIO_create_file_metadata
!EOC
!------------------------------------------------------------------------------
!BOP
! !ROUTINE: LIS_PFIO_create_routing_metadata
! \label{LIS_PFIO_create_routing_metadata}
! 
! !INTERFACE: LIS_PFIO_create_routing_metadata
   subroutine LIS_PFIO_create_routing_metadata(n, PFIOmodel_idx, group, outInterval, &
         nsoillayers, lyrthk, model_name)
! !USES: 
!
! !INPUTS PARAMETERS: 
      integer,          intent(in) :: n 
      integer,          intent(in) :: PFIOmodel_idx 
      integer,          intent(in) :: group
      real,             intent(in) :: outInterval
      integer,          intent(in) :: nsoillayers
      real,             intent(in) :: lyrthk(nsoillayers)
      character(len=*), intent(in) :: model_name
! 
! !DESCRIPTION: 
!  This routine creates a PFIO object that contains metadata to be used
!  in a netCDF file.
!  The arguments are: 
!  \begin{description}
!    \item[n] index of the nest
!    \item[ftn] file unit for the output file
!    \item[ftn\_stats] file unit for the output statistics file
!    \item[outInterval]   history output frequency
!    \item[nsoillayers]  Number of soil layers
!    \item[lyrthk]   Thickness of soil layers
!    \item[model\_name] Name of the model that generates the output
!  \end{description}
!
!   The routines invoked are: 
!   \begin{description}
!   \item[defineNETCDFheadervar](\ref{PFIO_define_variable_header})
!     writes the required headers for a single variable
!   \item[LIS\_verify](\ref{LIS_verify})
!     call to check if the return value is valid or not.
!   \end{description}
!
! !LOCAL VARIABLES:
      Type(Variable)          :: v

      integer                 :: dimID(4)
      integer                 :: tdimID,xtimeID,ensID
      integer                 :: t,c,r,i,index1, m
      real, allocatable       :: ensval(:) 
      character(len=8)        :: xtime_begin_date
      character(len=6)        :: xtime_begin_time
      character(len=50)       :: xtime_units
      character(len=50)       :: xtime_timeInc
      integer                 :: iret
      integer                 :: group_
      character(len=100)      :: model_name_
      integer                 :: status, gindex
      ! Note that the fix to add lat/lon to the NETCDF output will output
      ! undefined values for the water points. 
      character(len=8)        :: date
      character(len=10)       :: time
      character(len=90)       :: name_dims(4)
      character(len=5)        :: zone
      integer, dimension(8)   :: values
      type(LIS_metadataEntry), pointer :: dataEntry
      real                    :: time_data(1)
      REAL                    :: SOUTH_WEST_CORNER_LAT
      REAL                    :: SOUTH_WEST_CORNER_LON
      integer                 :: ierr
      integer                 :: vcol_id
!EOP
!---------------------------------------------------------------------------------------------
!BOC
      !--> ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      !--> We need to STOP if the "1d tilespace" configuration is selected.
      !--> PFIO cannot handle the way local tiles are mapped into the global one.
      !--> ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      IF (LIS_rc%wopt.eq."1d tilespace") THEN
         ierr = 1
         call LIS_verify(ierr,'The 1d tilespace configuration is not supported with PFIO!')
      ENDIF

      group_ = group
      model_name_ = model_name

      call LIS_rescaleCount(n, group_)

      call date_and_time(date,time,zone,values)

      pfio_vmin          = 0.0 ! -MAPL_UNDEF
      pfio_vmax          = 0.0 !  MAPL_UNDEF
      pfio_missing_value = LIS_rc%udef ! MAPL_UNDEF
      pfio_fill_value    = LIS_rc%udef ! MAPL_UNDEF
      pfio_valid_range   = (/-MAPL_UNDEF, MAPL_UNDEF/)

      !----------------------------
      ! Longitude and Latitude data
      !----------------------------

      allocate(xlat)
      allocate(xlong)

      xlat%short_name = "lat"
      xlat%long_name = "latitude"
      xlat%standard_name = "latitude"
      xlat%units = "degree_north"
      xlat%nunits = 1
      xlat%format = 'F'
      xlat%form = 1
      xlat%vlevels = 1
      xlat%timeAvgOpt = 0
      xlat%selectOpt = 1
      xlat%minMaxOpt = 0
      xlat%stdOpt = 0
      allocate(xlat%modelOutput(1,LIS_rc%nroutinggrid(n)*LIS_rc%nensem(n), xlat%vlevels))
      allocate(xlat%count(1,xlat%vlevels))
      xlat%count = 1
      allocate(xlat%unittypes(1))
      xlat%unittypes(1) = "degree_north"
      xlat%valid_min = 0.0
      xlat%valid_max = 0.0

      xlong%short_name = "lon"
      xlong%long_name = "longitude"
      xlong%standard_name = "longitude"
      xlong%units = "degree_east"
      xlong%nunits = 1
      xlong%format = 'F'
      xlong%form = 1
      xlong%vlevels = 1
      xlong%timeAvgOpt = 0
      xlong%selectOpt = 1
      xlong%minMaxOpt = 0
      xlong%stdOpt = 0
      allocate(xlong%modelOutput(1,LIS_rc%nroutinggrid(n)*LIS_rc%nensem(n),xlong%vlevels))
      allocate(xlong%count(1,xlong%vlevels))
      xlong%count = 1
      allocate(xlong%unittypes(1))
      xlong%unittypes(1) = "degree_east"
      xlong%valid_min = 0.0
      xlong%valid_max = 0.0

      ! Write variable data (both non-model and model-based):
      do i=1,LIS_rc%nroutinggrid(n)
         do m=1,LIS_rc%nensem(n)
            t = m+(i-1)*LIS_rc%nensem(n)

            c = LIS_routing(n)%tile(t)%col
            r = LIS_routing(n)%tile(t)%row

            xlat%modelOutput(1,t,1)  = LIS_domain(n)%lat(c+(r-1)*LIS_rc%lnc(n))
            xlong%modelOutput(1,t,1) = LIS_domain(n)%lon(c+(r-1)*LIS_rc%lnc(n))
         enddo
      enddo

      ! Write variable data (both non-model and model-based):
      if (LIS_rc%wopt.eq."2d ensemble gridspace") then
         allocate(ensval(LIS_rc%nensem(n)))
         do i = 1, LIS_rc%nensem(n)
            ensval(i) = float(i)
         end do
      endif

      ! defining time field
      write(xtime_units, 200) LIS_rc%yr, LIS_rc%mo, LIS_rc%da, LIS_rc%hr, LIS_rc%mn, LIS_rc%ss
      200   format ('minutes since ',I4.4,'-',I2.2,'-',I2.2,' ',I2.2,':', I2.2,':',I2.2)
      write(xtime_begin_date, fmt='(I4.4,I2.2,I2.2)') LIS_rc%yr, LIS_rc%mo, LIS_rc%da
      write(xtime_begin_time, fmt='(I2.2,I2.2,I2.2)') LIS_rc%hr, LIS_rc%mn, LIS_rc%ss
      write(xtime_timeInc, fmt='(I20)')  nint(outInterval)

      COL_LOOPS: DO vcol_id = 1, LIS_rc%n_vcollections
         write(LIS_logunit,'(a,i2,a1,i2)')'[INFO-PFIO] Create file metadata of collection ', vcol_id, '/', LIS_rc%n_vcollections
         name_dims(:) = ''

         !------------------
         ! Define dimensions
         !------------------
         if (LIS_rc%wopt.eq."1d tilespace") then 
            write(LIS_logunit,*) '[ERR] 1d tilespace output for routing models'
            write(LIS_logunit,*) '[ERR] is not supported currently'
            call LIS_endrun()
         elseif ( (LIS_rc%wopt.eq."2d gridspace") .OR. (LIS_rc%wopt.eq."2d ensemble gridspace") )then 
            name_dims(1) = 'east_west'
            name_dims(2) = 'north_south'
            call PFIO_bundle%fmd(n,vcol_id,PFIOmodel_idx)%add_dimension(TRIM(name_dims(1)), LIS_rc%gnc(n), rc=status) 
            call PFIO_bundle%fmd(n,vcol_id,PFIOmodel_idx)%add_dimension(TRIM(name_dims(2)), LIS_rc%gnr(n), rc=status)

            if (LIS_rc%wopt.eq."2d ensemble gridspace") then 
               name_dims(3) = 'ensemble'
               call PFIO_bundle%fmd(n,vcol_id,PFIOmodel_idx)%add_dimension(TRIM(name_dims(3)), LIS_rc%nensem(n), rc=status)

               v = Variable(type=PFIO_REAL32, dimensions='ensemble')
               call v%add_attribute("units", "ensemble number")
               call v%add_attribute("long_name", "Ensemble numbers")
               call v%add_const_value(UnlimitedEntity(ensval))
               call PFIO_bundle%fmd(n,vcol_id,PFIOmodel_idx)%add_variable('ensemble', v)
            endif
         endif

         ! LIS output is always writing output for a single time record
         call PFIO_bundle%fmd(n,vcol_id,PFIOmodel_idx)%add_dimension('time', pFIO_UNLIMITED, rc=status)

         !------------------------------------------
         ! Define variables and variables attributes
         !------------------------------------------

         ! Time: This is the initial step. 
         !       The attributes will be overwritten each time a file is created.
         time_data = 0.0
         v = Variable(type=PFIO_REAL32, dimensions='time')
         call v%add_attribute("units", trim(xtime_units))
         call v%add_attribute("long_name", "time")
         call v%add_attribute("time_increment", trim(adjustl(xtime_timeInc)))
         call v%add_attribute("begin_date", xtime_begin_date)
         call v%add_attribute("begin_time", xtime_begin_time)
         !call v%add_const_value(UnlimitedEntity( time_data ))
         call PFIO_bundle%fmd(n,vcol_id,PFIOmodel_idx)%add_variable('time', v)

         ! Pointer to header information
         dataEntry => LIS_histData(n)%head_routing_list

         total_num_2Dfields(PFIOmodel_idx) = 0
         total_num_3Dfields(PFIOmodel_idx) = 0
         total_num_4Dfields(PFIOmodel_idx) = 0

         !if ( (LIS_rc%wopt.eq."2d gridspace") .OR. (LIS_rc%wopt.eq."2d ensemble gridspace") )then
         call PFIO_define_variable_header(n, PFIOmodel_idx, vcol_id, xlat,  name_dims, non_model_fields = 1)
         call PFIO_define_variable_header(n, PFIOmodel_idx, vcol_id, xlong, name_dims, non_model_fields = 2)
         !endif

         do while ( associated(dataEntry) )
            call PFIO_define_variable_header(n, PFIOmodel_idx, vcol_id, dataEntry, name_dims)
            dataEntry => dataEntry%next
         enddo

         !------------------
         ! Global attributes
         !------------------
         call PFIO_bundle%fmd(n,vcol_id,PFIOmodel_idx)%add_attribute("missing_value", pfio_missing_value) ! LIS_rc%udef)
         call PFIO_bundle%fmd(n,vcol_id,PFIOmodel_idx)%add_attribute("NUM_SOIL_LAYERS", nsoillayers)
         call PFIO_bundle%fmd(n,vcol_id,PFIOmodel_idx)%add_attribute("SOIL_LAYER_THICKNESSES", lyrthk)
         call PFIO_bundle%fmd(n,vcol_id,PFIOmodel_idx)%add_attribute("title", "LIS land surface model output")
         call PFIO_bundle%fmd(n,vcol_id,PFIOmodel_idx)%add_attribute("institution", trim(LIS_rc%institution))
         call PFIO_bundle%fmd(n,vcol_id,PFIOmodel_idx)%add_attribute("source", trim(model_name_))
         call PFIO_bundle%fmd(n,vcol_id,PFIOmodel_idx)%add_attribute("history",  & 
            "created on date: "//date(1:4)//"-"//date(5:6)//"-"// date(7:8)//"T"//time(1:2)//":"//time(3:4)//":"//time(5:10))
         call PFIO_bundle%fmd(n,vcol_id,PFIOmodel_idx)%add_attribute("references", "Kumar_etal_EMS_2006, Peters-Lidard_etal_ISSE_2007")
         call PFIO_bundle%fmd(n,vcol_id,PFIOmodel_idx)%add_attribute("conventions", "CF-1.6")
         call PFIO_bundle%fmd(n,vcol_id,PFIOmodel_idx)%add_attribute("comment", "website: http://lis.gsfc.nasa.gov/")

         ! Grid information

         !============================================================================
         ! Make sure that all the processes have the same value the lower left corner.
         ! Broadcast the the root process values to all the other processes.
         ! This is only needed to write global attributes in the file.
         !============================================================================
         IF (LIS_masterproc) THEN
            SOUTH_WEST_CORNER_LAT = LIS_rc%gridDesc(n,4)
            SOUTH_WEST_CORNER_LON = LIS_rc%gridDesc(n,5)
         ELSE
            SOUTH_WEST_CORNER_LAT = -9999.0
            SOUTH_WEST_CORNER_LON = -9999.0
         ENDIF
         call MPI_Bcast(SOUTH_WEST_CORNER_LAT,1, MPI_REAL, 0, LIS_mpi_comm, ierr)
         call MPI_Bcast(SOUTH_WEST_CORNER_LON,1, MPI_REAL, 0, LIS_mpi_comm, ierr)

         if (LIS_rc%lis_map_proj.eq."latlon") then   ! latlon
            call PFIO_bundle%fmd(n,vcol_id,PFIOmodel_idx)%add_attribute("MAP_PROJECTION", "EQUIDISTANT CYLINDRICAL")
            call PFIO_bundle%fmd(n,vcol_id,PFIOmodel_idx)%add_attribute("SOUTH_WEST_CORNER_LAT", SOUTH_WEST_CORNER_LAT)
            call PFIO_bundle%fmd(n,vcol_id,PFIOmodel_idx)%add_attribute("SOUTH_WEST_CORNER_LON", SOUTH_WEST_CORNER_LON)
            call PFIO_bundle%fmd(n,vcol_id,PFIOmodel_idx)%add_attribute("DX", LIS_rc%gridDesc(n,9))
            call PFIO_bundle%fmd(n,vcol_id,PFIOmodel_idx)%add_attribute("DY", LIS_rc%gridDesc(n,10))       
         elseif (LIS_rc%lis_map_proj.eq."mercator") then 
            call PFIO_bundle%fmd(n,vcol_id,PFIOmodel_idx)%add_attribute("MAP_PROJECTION", "MERCATOR")
            call PFIO_bundle%fmd(n,vcol_id,PFIOmodel_idx)%add_attribute("SOUTH_WEST_CORNER_LAT", SOUTH_WEST_CORNER_LAT)
            call PFIO_bundle%fmd(n,vcol_id,PFIOmodel_idx)%add_attribute("SOUTH_WEST_CORNER_LON", SOUTH_WEST_CORNER_LON)
            call PFIO_bundle%fmd(n,vcol_id,PFIOmodel_idx)%add_attribute("TRUELAT1", LIS_rc%gridDesc(n,10))
            call PFIO_bundle%fmd(n,vcol_id,PFIOmodel_idx)%add_attribute("STANDARD_LON", LIS_rc%gridDesc(n,11))
            call PFIO_bundle%fmd(n,vcol_id,PFIOmodel_idx)%add_attribute("DX", LIS_rc%gridDesc(n,8))
            call PFIO_bundle%fmd(n,vcol_id,PFIOmodel_idx)%add_attribute("DY", LIS_rc%gridDesc(n,9))
         elseif (LIS_rc%lis_map_proj.eq."lambert") then ! lambert conformal
            call PFIO_bundle%fmd(n,vcol_id,PFIOmodel_idx)%add_attribute("MAP_PROJECTION",  "LAMBERT CONFORMAL")
            call PFIO_bundle%fmd(n,vcol_id,PFIOmodel_idx)%add_attribute("SOUTH_WEST_CORNER_LAT", SOUTH_WEST_CORNER_LAT)
            call PFIO_bundle%fmd(n,vcol_id,PFIOmodel_idx)%add_attribute("SOUTH_WEST_CORNER_LON", SOUTH_WEST_CORNER_LON)
            call PFIO_bundle%fmd(n,vcol_id,PFIOmodel_idx)%add_attribute("TRUELAT1", LIS_rc%gridDesc(n,10))
            call PFIO_bundle%fmd(n,vcol_id,PFIOmodel_idx)%add_attribute("TRUELAT2", LIS_rc%gridDesc(n,7))
            call PFIO_bundle%fmd(n,vcol_id,PFIOmodel_idx)%add_attribute("STANDARD_LON", LIS_rc%gridDesc(n,11))
            call PFIO_bundle%fmd(n,vcol_id,PFIOmodel_idx)%add_attribute("DX", LIS_rc%gridDesc(n,8))
            call PFIO_bundle%fmd(n,vcol_id,PFIOmodel_idx)%add_attribute("DY", LIS_rc%gridDesc(n,9))
         elseif (LIS_rc%lis_map_proj.eq."polar") then ! polar stereographic
            call PFIO_bundle%fmd(n,vcol_id,PFIOmodel_idx)%add_attribute("MAP_PROJECTION", "POLAR STEREOGRAPHIC")
            call PFIO_bundle%fmd(n,vcol_id,PFIOmodel_idx)%add_attribute("SOUTH_WEST_CORNER_LAT", SOUTH_WEST_CORNER_LAT)
            call PFIO_bundle%fmd(n,vcol_id,PFIOmodel_idx)%add_attribute("SOUTH_WEST_CORNER_LON", SOUTH_WEST_CORNER_LON)
            call PFIO_bundle%fmd(n,vcol_id,PFIOmodel_idx)%add_attribute("TRUELAT1", LIS_rc%gridDesc(n,10))
            call PFIO_bundle%fmd(n,vcol_id,PFIOmodel_idx)%add_attribute("ORIENT", LIS_rc%gridDesc(n,7))
            call PFIO_bundle%fmd(n,vcol_id,PFIOmodel_idx)%add_attribute("STANDARD_LON", LIS_rc%gridDesc(n,11))
            call PFIO_bundle%fmd(n,vcol_id,PFIOmodel_idx)%add_attribute("DX", LIS_rc%gridDesc(n,8))
            call PFIO_bundle%fmd(n,vcol_id,PFIOmodel_idx)%add_attribute("DY", LIS_rc%gridDesc(n,9))
         endif       

         ! Create the history collection
         PFIO_bundle%hist_id(n,vcol_id,PFIOmodel_idx) = o_Clients%add_hist_collection(PFIO_bundle%fmd(n,vcol_id,PFIOmodel_idx))
      ENDDO COL_LOOPS

      ! ----> Done with the creation of the PFIO object

      if ( LIS_rc%wopt.eq."2d ensemble gridspace" ) then
         deallocate(ensval)
      endif

   end subroutine LIS_PFIO_create_routing_metadata
!EOC
!------------------------------------------------------------------------------
!BOP
! !ROUTINE: PFIO_define_variable_header
! \label{PFIO_define_variable_header}
! 
! !INTERFACE: 
!
   subroutine PFIO_define_variable_header(n, PFIOmodel_idx, vcol_id, dataEntry, name_dims, non_model_fields)
!
! !USES: 

! !ARGUMENTS:     
      integer, intent(in)               :: n
      integer, intent(in)               :: PFIOmodel_idx
      integer, intent(in)               :: vcol_id
      type(LIS_metadataEntry), pointer  :: dataEntry
      character(len=*)                  :: name_dims(4)
      integer,   optional               :: non_model_fields
! 
! !DESCRIPTION: 
!    This routine writes the required NETCDF header for a single variable
! 
!   The arguments are: 
!   \begin{description}
!   \item[n]
!    index of the nest
!   \item[fmd]
!    PFIO file metadata
!   \item[dataEntry]
!    object containing the values and attributes of the variable to be 
!    written
!   \end{description}
!
!   The routines invoked are: 
!   \begin{description}
!   \item[LIS\_endrun](\ref{LIS_endrun})
!     call to abort program when a fatal error is detected. 
!   \item[LIS\_verify](\ref{LIS_verify})
!     call to check if the return value is valid or not.
!   \end{description}
!
! !LOCAL VARIABLES:
      integer            :: data_index
      integer            :: status, ierr
      integer            :: deflate_level
      character(len=100) :: short_name
      integer            :: fill_value

      character(len=100) :: vname_def, vname_opt1, vname_min, vname_max
      character(len=256) :: dim_names
      Type(Variable)     :: v_def, v_opt1, v_min, v_max
      REAL               :: vmin 
      REAL               :: vmax
      integer            :: nmodel_status
!EOP
!------------------------------------------------------------------------------
!BOC

      nmodel_status = 0
      if (present(non_model_fields)) then
         nmodel_status = non_model_fields
      endif

      data_index = dataEntry%index

      deflate_level = LIS_rc%nc_deflate_lvl

      if (dataEntry%selectOpt.eq.1) then 
         if (LIS_rc%wopt.eq."1d tilespace") then              
            if (dataEntry%vlevels.gt.1) then 
               name_dims(2) = trim(dataEntry%short_name)//'_profiles'
               call PFIO_bundle%fmd(n,vcol_id,PFIOmodel_idx)%add_dimension(trim(name_dims(2)), dataEntry%vlevels, rc=status)
            endif
         elseif(LIS_rc%wopt.eq."2d gridspace") then
            if(dataEntry%vlevels.gt.1) then
               name_dims(3) = trim(dataEntry%short_name)//'_profiles'
               call PFIO_bundle%fmd(n,vcol_id,PFIOmodel_idx)%add_dimension(trim(name_dims(3)), dataEntry%vlevels, rc=status)
            endif
         elseif(LIS_rc%wopt.eq."2d ensemble gridspace") then
            if(dataEntry%vlevels.gt.1) then
               name_dims(4) = trim(dataEntry%short_name)//'_profiles'
               call PFIO_bundle%fmd(n,vcol_id,PFIOmodel_idx)%add_dimension(trim(name_dims(4)), dataEntry%vlevels, rc=status)
            endif
         endif

         if (LIS_rc%wopt.eq."1d tilespace") then

            if (dataEntry%timeAvgOpt.eq.2) then 
               if (dataEntry%vlevels.gt.1) then 
                  dim_names = TRIM(name_dims(1))//','//TRIM(name_dims(2))
               else
                  dim_names = TRIM(name_dims(1))
               endif

               total_num_2Dfields(PFIOmodel_idx) = total_num_2Dfields(PFIOmodel_idx) + 2

               vname_def  = TRIM(dataEntry%short_name)//'_tavg'
               v_def  = Variable(type=PFIO_REAL32, dimensions=TRIM(dim_names), deflation=deflate_level)

               vname_opt1 = TRIM(dataEntry%short_name)//'_inst'
               v_opt1 = Variable(type=PFIO_REAL32, dimensions=TRIM(dim_names), deflation=deflate_level)

               if (dataEntry%minMaxOpt.gt.0) then 
                  total_num_2Dfields(PFIOmodel_idx) = total_num_2Dfields(PFIOmodel_idx) + 2
                  vname_min = trim(dataEntry%short_name)//"_min"
                  v_min = Variable(type=PFIO_REAL32, dimensions=TRIM(dim_names), deflation=deflate_level)

                  vname_max = trim(dataEntry%short_name)//"_max"
                  v_max = Variable(type=PFIO_REAL32, dimensions=TRIM(dim_names), deflation=deflate_level)
               endif
            else
               if (nmodel_status == 0) then
                  if (dataEntry%timeAvgOpt.eq.0) then
                     short_name = trim(dataEntry%short_name)//'_inst'
                  elseif (dataEntry%timeAvgOpt.eq.1) then
                     short_name = trim(dataEntry%short_name)//'_tavg'
                  elseif (dataEntry%timeAvgOpt.eq.3) then
                     short_name = trim(dataEntry%short_name)//'_acc'
                  endif
               else
                  short_name = trim(dataEntry%short_name)
               endif

               if (dataEntry%vlevels.gt.1) then 
                  dim_names = TRIM(name_dims(1))//','//TRIM(name_dims(2))
               else
                  dim_names = TRIM(name_dims(1))
               endif

               total_num_2Dfields(PFIOmodel_idx) = total_num_2Dfields(PFIOmodel_idx) + 1
               vname_def = trim(short_name)
               v_def     = Variable(type=PFIO_REAL32, dimensions=TRIM(dim_names), deflation=deflate_level)

               if (dataEntry%minMaxOpt.gt.0) then 
                  total_num_2Dfields(PFIOmodel_idx) = total_num_2Dfields(PFIOmodel_idx) + 2
                  dim_names = TRIM(name_dims(1))
                  vname_min = trim(short_name)//"_min"
                  vname_max = trim(short_name)//"_max"
                  v_min     = Variable(type=PFIO_REAL32, dimensions=TRIM(dim_names), deflation=deflate_level)
                  v_max     = Variable(type=PFIO_REAL32, dimensions=TRIM(dim_names), deflation=deflate_level)
               endif
            endif
         elseif (LIS_rc%wopt.eq."2d gridspace") then 
            if (dataEntry%timeAvgOpt.eq.2) then 
               if (dataEntry%vlevels.gt.1) then 
                  dim_names  = TRIM(name_dims(1))//','//TRIM(name_dims(2))//','//TRIM(name_dims(3))
               else
                  dim_names = TRIM(name_dims(1))//','//TRIM(name_dims(2))
               endif

               total_num_3Dfields(PFIOmodel_idx) = total_num_3Dfields(PFIOmodel_idx) + 2

               vname_def  = TRIM(dataEntry%short_name)//'_tavg'
               v_def      = Variable(type=PFIO_REAL32, dimensions=TRIM(dim_names), deflation=deflate_level)

               vname_opt1 = TRIM(dataEntry%short_name)//'_inst'
               v_opt1     = Variable(type=PFIO_REAL32, dimensions=TRIM(dim_names), deflation=deflate_level)

               if (dataEntry%minMaxOpt.gt.0) then 
                  total_num_3Dfields(PFIOmodel_idx) = total_num_3Dfields(PFIOmodel_idx) + 2
                  vname_min = trim(dataEntry%short_name)//"_min"
                  vname_max = trim(dataEntry%short_name)//"_max"
                  v_min     = Variable(type=PFIO_REAL32, dimensions=TRIM(dim_names), deflation=deflate_level)
                  v_max     = Variable(type=PFIO_REAL32, dimensions=TRIM(dim_names), deflation=deflate_level)
               endif
            else
               if (nmodel_status == 0) then
                  if (dataEntry%timeAvgOpt.eq.0) then
                     short_name = trim(dataEntry%short_name)//'_inst'
                  elseif (dataEntry%timeAvgOpt.eq.1) then
                     short_name = trim(dataEntry%short_name)//'_tavg'
                  elseif (dataEntry%timeAvgOpt.eq.3) then
                     short_name = trim(dataEntry%short_name)//'_acc'
                  endif
               else
                  short_name = trim(dataEntry%short_name)
               endif

               if (dataEntry%vlevels.gt.1) then 
                  dim_names = TRIM(name_dims(1))//','//TRIM(name_dims(2))//','//TRIM(name_dims(3))
               else
                  ! lat/lon fields will write in 1D  
                  if (LIS_rc%nlatlon_dimensions == '1D') then
                     if (nmodel_status.eq.1) then
                        dim_names = TRIM(name_dims(2))
                     elseif (nmodel_status.eq.2) then
                        dim_names = TRIM(name_dims(1))
                     else
                        dim_names = TRIM(name_dims(1))//','//TRIM(name_dims(2))
                     endif
                  else
                     dim_names = TRIM(name_dims(1))//','//TRIM(name_dims(2))
                  endif
               endif

               total_num_3Dfields(PFIOmodel_idx) = total_num_3Dfields(PFIOmodel_idx) + 1
               vname_def = trim(short_name)
               v_def = Variable(type=PFIO_REAL32, dimensions=TRIM(dim_names), deflation=deflate_level)

               if (dataEntry%minMaxOpt.gt.0) then 
                  total_num_3Dfields(PFIOmodel_idx) = total_num_3Dfields(PFIOmodel_idx) + 2
                  vname_min = trim(short_name)//"_min"
                  vname_max = trim(short_name)//"_max"

                  v_min = Variable(type=PFIO_REAL32, dimensions=TRIM(dim_names), deflation=deflate_level)
                  v_max = Variable(type=PFIO_REAL32, dimensions=TRIM(dim_names), deflation=deflate_level)

               endif
            endif
         elseif (LIS_rc%wopt.eq."2d ensemble gridspace") then 

            if (dataEntry%timeAvgOpt.eq.2) then 

               if (dataEntry%vlevels.gt.1) then 
                  dim_names  = TRIM(name_dims(1))//','//TRIM(name_dims(2))//','//TRIM(name_dims(3))//','//TRIM(name_dims(4))

                  total_num_4Dfields(PFIOmodel_idx) = total_num_4Dfields(PFIOmodel_idx) + 2
                  vname_def  = TRIM(dataEntry%short_name)//'_tavg'
                  vname_opt1 = TRIM(dataEntry%short_name)//'_inst'
                  v_def      = Variable(type=PFIO_REAL32, dimensions=TRIM(dim_names), deflation=deflate_level)
                  v_opt1     = Variable(type=PFIO_REAL32, dimensions=TRIM(dim_names), deflation=deflate_level)

                  if (dataEntry%minMaxOpt.gt.0) then 
                     total_num_4Dfields(PFIOmodel_idx) = total_num_4Dfields(PFIOmodel_idx) + 2
                     vname_min = trim(dataEntry%short_name)//"_min"
                     vname_max = trim(dataEntry%short_name)//"_max"
                     v_min     = Variable(type=PFIO_REAL32, dimensions=TRIM(dim_names), deflation=deflate_level)
                     v_max     = Variable(type=PFIO_REAL32, dimensions=TRIM(dim_names), deflation=deflate_level)
                  endif
               else
                  short_name = trim(dataEntry%short_name)
                  if (nmodel_status == 0) then
                     dim_names  = TRIM(name_dims(1))//','//TRIM(name_dims(2))//','//TRIM(name_dims(3))

                     total_num_4Dfields(PFIOmodel_idx) = total_num_4Dfields(PFIOmodel_idx) + 2
                     vname_def  = TRIM(dataEntry%short_name)//'_tavg'
                     vname_opt1 = TRIM(dataEntry%short_name)//'_inst'
                     v_def      = Variable(type=PFIO_REAL32, dimensions=TRIM(dim_names), deflation=deflate_level)
                     v_opt1     = Variable(type=PFIO_REAL32, dimensions=TRIM(dim_names), deflation=deflate_level)

                  else
                     vname_def = trim(short_name)
                     ! lat/lon fields will write in 1D  
                     if (LIS_rc%nlatlon_dimensions == '1D') then
                        if (nmodel_status.eq.1) then
                           dim_names = TRIM(name_dims(2))
                        elseif (nmodel_status.eq.2) then
                           dim_names = TRIM(name_dims(1))
                        else
                           dim_names = TRIM(name_dims(1))//','//TRIM(name_dims(2))
                        endif
                     else
                        dim_names = TRIM(name_dims(1))//','//TRIM(name_dims(2))
                     endif
                     v_def = Variable(type=PFIO_REAL32, dimensions=TRIM(dim_names), deflation=deflate_level)
                  endif
                  if (dataEntry%minMaxOpt.gt.0) then 
                     total_num_3Dfields(PFIOmodel_idx) = total_num_3Dfields(PFIOmodel_idx) + 2
                     vname_min = trim(dataEntry%short_name)//"_min"
                     vname_max = trim(dataEntry%short_name)//"_max"
                     v_min     = Variable(type=PFIO_REAL32, dimensions=TRIM(dim_names), deflation=deflate_level)
                     v_max     = Variable(type=PFIO_REAL32, dimensions=TRIM(dim_names), deflation=deflate_level)
                  endif
               endif
            else
               if (nmodel_status > 0) then
                  short_name = trim(dataEntry%short_name)
               else
                  if (dataEntry%timeAvgOpt.eq.0) then
                     short_name = trim(dataEntry%short_name)//'_inst'
                  elseif (dataEntry%timeAvgOpt.eq.1) then
                     short_name = trim(dataEntry%short_name)//'_tavg'
                  elseif (dataEntry%timeAvgOpt.eq.3) then
                     short_name = trim(dataEntry%short_name)//'_acc'
                  endif
               endif

               vname_def = trim(short_name)

               if (dataEntry%vlevels.gt.1) then 
                  total_num_4Dfields(PFIOmodel_idx) = total_num_4Dfields(PFIOmodel_idx) + 1
                  dim_names = TRIM(name_dims(1))//','//TRIM(name_dims(2))//','//TRIM(name_dims(3))//','//TRIM(name_dims(4))
                  v_def = Variable(type=PFIO_REAL32, dimensions=TRIM(dim_names), deflation=deflate_level)

               else
                  if (nmodel_status == 0) then
                     total_num_4Dfields(PFIOmodel_idx) = total_num_4Dfields(PFIOmodel_idx) + 1
                     dim_names = TRIM(name_dims(1))//','//TRIM(name_dims(2))//','//TRIM(name_dims(3))
                  else
                     vname_def = trim(short_name)
                     ! lat/lon fields will write in 1D  
                     if (LIS_rc%nlatlon_dimensions == '1D') then
                        if (nmodel_status.eq.1) then
                           dim_names = TRIM(name_dims(2))
                        elseif (nmodel_status.eq.2) then
                           dim_names = TRIM(name_dims(1))
                        endif
                     else
                        dim_names = TRIM(name_dims(1))//','//TRIM(name_dims(2))
                     endif
                  end if
                  v_def = Variable(type=PFIO_REAL32, dimensions=TRIM(dim_names), deflation=deflate_level)

                  if (dataEntry%minMaxOpt.gt.0) then 
                     total_num_4Dfields(PFIOmodel_idx) = total_num_4Dfields(PFIOmodel_idx) + 2
                     vname_min = trim(short_name)//"_min"
                     vname_max = trim(short_name)//"_max"
                     v_min     = Variable(type=PFIO_REAL32, dimensions=TRIM(dim_names), deflation=deflate_level)
                     v_max     = Variable(type=PFIO_REAL32, dimensions=TRIM(dim_names), deflation=deflate_level)
                  endif
               endif
            endif
         endif

         call v_def%add_attribute("units",  trim(dataEntry%units))
         call v_def%add_attribute("standard_name", trim(dataEntry%standard_name))
         call v_def%add_attribute("long_name", trim(dataEntry%long_name))
         call v_def%add_attribute("scale_factor", 1.0)
         call v_def%add_attribute("add_offset", 0.0)
         call v_def%add_attribute("missing_value", pfio_missing_value) 
         call v_def%add_attribute("_FillValue", pfio_fill_value) 
         !call v_def%add_attribute('valid_range', pfio_valid_range)
         call v_def%add_attribute("vmin", pfio_vmin)
         call v_def%add_attribute("vmax", pfio_vmax)
         call PFIO_bundle%fmd(n,vcol_id,PFIOmodel_idx)%add_variable(TRIM(vname_def), v_def)

         if (dataEntry%timeAvgOpt.eq.2) then
            call v_opt1%add_attribute("units", trim(dataEntry%units))
            call v_opt1%add_attribute("standard_name", trim(dataEntry%standard_name))
            call v_opt1%add_attribute("long_name", trim(dataEntry%long_name))
            call v_opt1%add_attribute("scale_factor", 1.0)
            call v_opt1%add_attribute("add_offset", 0.0)
            call v_opt1%add_attribute("missing_value", pfio_missing_value) 
            call v_opt1%add_attribute("_FillValue", pfio_fill_value) 
            !call v_opt1%add_attribute('valid_range', pfio_valid_range)
            call v_opt1%add_attribute("vmin", pfio_vmin)
            call v_opt1%add_attribute("vmax", pfio_vmax)
            call PFIO_bundle%fmd(n,vcol_id,PFIOmodel_idx)%add_variable(TRIM(vname_opt1), v_opt1)
         endif

         if (dataEntry%minMaxOpt.gt.0) then
            ! Min metadata
            call v_min%add_attribute("units", trim(dataEntry%units))
            call v_min%add_attribute("standard_name", trim(dataEntry%standard_name))
            call v_min%add_attribute("long_name", trim(dataEntry%long_name))
            call v_min%add_attribute("scale_factor", 1.0)
            call v_min%add_attribute("add_offset", 0.0)
            call v_min%add_attribute("missing_value", pfio_missing_value) 
            call v_min%add_attribute("_FillValue", pfio_fill_value) 
            !call v_min%add_attribute('valid_range', pfio_valid_range)
            call v_min%add_attribute("vmin", pfio_vmin)
            call v_min%add_attribute("vmax", pfio_vmax)
            call PFIO_bundle%fmd(n,vcol_id,PFIOmodel_idx)%add_variable(TRIM(vname_min), v_min)

            ! Max metadata
            call v_max%add_attribute("units", trim(dataEntry%units))
            call v_max%add_attribute("standard_name", trim(dataEntry%standard_name))
            call v_max%add_attribute("long_name", trim(dataEntry%long_name))
            call v_max%add_attribute("scale_factor", 1.0)
            call v_max%add_attribute("add_offset", 0.0)
            call v_max%add_attribute("missing_value", pfio_missing_value) 
            call v_max%add_attribute("_FillValue", pfio_fill_value) 
            !call v_max%add_attribute('valid_range', pfio_valid_range)
            call v_max%add_attribute("vmin", pfio_vmin)
            call v_max%add_attribute("vmax", pfio_vmax)
            call PFIO_bundle%fmd(n,vcol_id,PFIOmodel_idx)%add_variable(TRIM(vname_max), v_max)
         endif
      endif

   end subroutine PFIO_define_variable_header
!EOC
!------------------------------------------------------------------------------
!BOP
! !ROUTINE: LIS_PFIO_write_data
! \label{LIS_PFIO_write_data}
!
! !INTERFACE: 
   subroutine LIS_PFIO_write_data(n, PFIOmodel_idx, vcol_id, file_name, outInterval, group)

      integer,          intent(in) :: n 
      integer,          intent(in) :: PFIOmodel_idx 
      integer,          intent(in) :: vcol_id 
      character(len=*), intent(in) :: file_name
      real,             intent(in) :: outInterval
      integer,   optional     :: group
! 
! !DESCRIPTION: 
!  This routine writes variables to a NETCDF file
!
! !LOCAL VARIABLES:
      real                    :: time_data(1)
      integer                 :: i,k,m,t, group_, status, rc
      character(len=8)        :: xtime_begin_date
      character(len=6)        :: xtime_begin_time
      character(len=50)       :: xtime_units
      character(len=50)       :: xtime_timeInc
      character(len=8)        :: date
      character(len=10)       :: time
      character(len=5)        :: zone
      integer, dimension(8)   :: values
      Type(Variable)          :: v
      type(StringVariableMap) :: var_map
      type(LIS_metadataEntry), pointer :: dataEntry
      type(local_2Dfield_var), allocatable :: local_var2D(:)
      type(local_3Dfield_var), allocatable :: local_var3D(:)
      type(local_4Dfield_var), allocatable :: local_var4D(:)
      integer                 :: i1, i2, j1, j2
      type(ArrayReference)    :: ref
      integer                 :: global_dim(2)
      integer                 :: r, c, index1
!EOP
!---------------------------------------------------------------------------------------------
!BOC
      write(LIS_logunit,'(a,i2,a,a)')'[INFO-PFIO] Writing model output [',PFIOmodel_idx,'] to: ', TRIM(file_name)

      ! Update the time variable
      !-------------------------
      call date_and_time(date, time, zone, values)
      write(xtime_units, 200) LIS_rc%yr, LIS_rc%mo, LIS_rc%da, LIS_rc%hr, LIS_rc%mn, LIS_rc%ss
      200   format ('minutes since ',I4.4,'-',I2.2,'-',I2.2,' ',I2.2,':', I2.2,':',I2.2)
      write(xtime_begin_date, fmt='(I4.4,I2.2,I2.2)') LIS_rc%yr, LIS_rc%mo, LIS_rc%da
      write(xtime_begin_time, fmt='(I2.2,I2.2,I2.2)') LIS_rc%hr, LIS_rc%mn, LIS_rc%ss
      write(xtime_timeInc, fmt='(I20)')  nint(outInterval)

      time_data = (/ 1.0 /)
      v = Variable(type=PFIO_REAL32, dimensions='time')
      call v%add_attribute("units",trim(xtime_units))
      call v%add_attribute("long_name","time")
      call v%add_attribute("time_increment",trim(adjustl(xtime_timeInc)))
      call v%add_attribute("begin_date",xtime_begin_date)
      call v%add_attribute("begin_time",xtime_begin_time)
      !call v%add_const_value(UnlimitedEntity( time_data ))
      call PFIO_bundle%fmd(n,vcol_id,PFIOmodel_idx)%add_variable('time', v)
      call var_map%insert('time', v)
      call o_Clients%modify_metadata(PFIO_bundle%hist_id(n,vcol_id,PFIOmodel_idx), var_map=var_map, rc=status)

      call o_Clients%set_optimal_server(nwriting=1)

      group_ = 1
      if(PRESENT(group)) group_ = group

      if (LIS_rc%wopt.eq."1d tilespace") then
         CALL get_LIS_interior_tile(n, i1, i2)
      elseif ( (LIS_rc%wopt.eq."2d gridspace") .OR. (LIS_rc%wopt.eq."2d ensemble gridspace") )then 
         CALL get_LIS_interior_grid(n, i1, i2, j1, j2)
      endif
      CALL get_LIS_globaldim(n, global_dim)


      SELECT CASE (group_)
      CASE(1) ! LSM output
         dataEntry => LIS_histData(n)%head_lsm_list
      CASE(2) ! ROUTING
         dataEntry => LIS_histData(n)%head_routing_list
      CASE(3) ! RTM
         dataEntry => LIS_histData(n)%head_rtm_list
      CASE(4) ! Irrigation
         dataEntry => LIS_histData(n)%head_irrig_list
      END SELECT

      SELECT CASE(LIS_rc%wopt)
      CASE("1d tilespace")
         ALLOCATE(local_var2D(total_num_2Dfields(PFIOmodel_idx)))
      CASE("2d gridspace")
         ALLOCATE(local_var3D(total_num_3Dfields(PFIOmodel_idx)))
      CASE("2d ensemble gridspace")
         ALLOCATE(local_var3D(total_num_3Dfields(PFIOmodel_idx)))
         ALLOCATE(local_var4D(total_num_4Dfields(PFIOmodel_idx)))
      END SELECT
      call allocate_local_var(n, PFIOmodel_idx, local_var2D, local_var3D, local_var4D, i1, i2, j1, j2)
      idx_field2d = 0
      idx_field3d = 0
      idx_field4d = 0

      ! Writing latitude/longitude
      !---------------------------
      !if ( (LIS_rc%wopt.eq."2d gridspace") .OR. (LIS_rc%wopt.eq."2d ensemble gridspace") )then 
      call PFIO_write_variable(n, PFIOmodel_idx, vcol_id, xlat,  &
         file_name, local_var2D, local_var3D, local_var4D, &
         i1, i2, j1, j2, non_model_fields=1)
      call PFIO_write_variable(n, PFIOmodel_idx, vcol_id, xlong, &
         file_name, local_var2D, local_var3D, local_var4D, &
         i1, i2, j1, j2, non_model_fields=2)
      !endif

      ! Writing the model fields
      !-------------------------
      do while ( associated(dataEntry) )
         call PFIO_write_variable(n, PFIOmodel_idx, vcol_id, dataEntry, file_name, &
            local_var2D, local_var3D, local_var4D, i1, i2, j1, j2)
         dataEntry => dataEntry%next
      enddo

      ! write in the file and close it
      call o_Clients%done_collective_stage()

      call o_Clients%post_wait()

      ! After writing reset the variables
      SELECT CASE(LIS_rc%wopt)
      CASE("1d tilespace")
         DO i= 1, total_num_2Dfields(PFIOmodel_idx)
            deallocate(local_var2D(i)%var2d)
         ENDDO
         deallocate(local_var2D)
      CASE("2d gridspace")
         DO i= 1, total_num_3Dfields(PFIOmodel_idx)
            deallocate(local_var3D(i)%var3d)
         ENDDO
         deallocate(local_var3D)
      CASE("2d ensemble gridspace")
         DO i= 1, total_num_3Dfields(PFIOmodel_idx)
            deallocate(local_var3D(i)%var3d)
         ENDDO
         deallocate(local_var3D)
         DO i= 1, total_num_4Dfields(PFIOmodel_idx)
            deallocate(local_var4D(i)%var4d)
         ENDDO
         deallocate(local_var4D)
      END SELECT
      !CALL deallocate_local_var(PFIOmodel_idx, local_var2D, local_var3D, local_var4D)
      !IF (ALLOCATED(local_var2D)) DEALLOCATE(local_var2D)
      !IF (ALLOCATED(local_var3D)) DEALLOCATE(local_var3D)
      !IF (ALLOCATED(local_var4D)) DEALLOCATE(local_var4D)

      call LIS_resetOutputVars(n, group_)

   end subroutine LIS_PFIO_write_data
!EOC
!------------------------------------------------------------------------------
!BOP
! !ROUTINE: PFIO_write_variable
! \label{PFIO_write_variable}
!
! !INTERFACE: 
   subroutine PFIO_write_variable(n, PFIOmodel_idx, vcol_id, dataEntry, file_name, &
         local_var2D, local_var3D, local_var4D, &
         i1, i2, j1, j2, non_model_fields)

      integer,   intent(in)   :: n 
      integer,   intent(in)   :: PFIOmodel_idx 
      integer,   intent(in)   :: vcol_id 
      integer,   intent(in)   :: i1, i2, j1, j2
      type(local_2Dfield_var), intent(inOut) :: local_var2D(:)
      type(local_3Dfield_var), intent(inOut) :: local_var3D(:)
      type(local_4Dfield_var), intent(inOut) :: local_var4D(:)
      type(LIS_metadataEntry), pointer :: dataEntry
      character(len=*), intent(in) :: file_name
      integer, optional, intent(in) :: non_model_fields
! 
! !DESCRIPTION: 
!  This routine writes a single variable to a NETCDF file
!  The arguments are: 
!  \begin{description}
!    \item[n] index of the nest
!   \item[dataEntry]
!    object containing the values and attributes of the variable to be 
!    written
!  \end{description}
!
!   The routines invoked are: 
!   \begin{description}
!   \item[LIS\_writevar\_netcdf](\ref{LIS_writevar_netcdf})
!     writes a variable into a netcdf formatted file. 
!   \end{description}
!
! !LOCAL VARIABLES:
      integer       :: i,k,m,t, nlev
      character(len=90) :: var_name
      character(len=90) :: var_name2
      integer       :: nmodel_status
!EOP    
!------------------------------------------------------------------------------
!BOC

      nmodel_status = 0
      if (present(non_model_fields)) then
         nmodel_status = non_model_fields
      endif

      IF_selectOpt: if (dataEntry%selectOpt.eq.1) then

         m = 1
         IF (nmodel_status == 0) THEN
            do t=1,LIS_rc%ntiles(n)
               m = LIS_domain(n)%tile(t)%sftype
               do k=1,dataEntry%vlevels
                  if (dataEntry%count(t,k).gt.0) then 
                     SELECT CASE(dataEntry%timeAvgOpt)
                     CASE(3)
                        continue   
                     CASE (1, 2)
                        dataEntry%modelOutput(1,t,k) = dataEntry%modelOutput(1,t,k)/ dataEntry%count(t,k)
                     CASE DEFAULT
                        continue   
                     END SELECT
                  else
                     dataEntry%modelOutput(1,t,k) = pfio_missing_value
                  endif
               enddo
            enddo
         ENDIF

         if (nmodel_status == 0) then
            SELECT CASE(dataEntry%timeAvgOpt)
            CASE(0)
               var_name = trim(dataEntry%short_name)//'_inst'
            CASE(1)
               var_name = trim(dataEntry%short_name)//'_tavg'
            CASE(3)
               var_name = trim(dataEntry%short_name)//'_acc'
            END SELECT
         else
            var_name = trim(dataEntry%short_name)
         endif

         ! accumulated values
         ! time-averaged values and instantaneous values
         nlev = dataEntry%vlevels
         if (dataEntry%timeAvgOpt.eq.2) then 
            CALL increment_field_counter()
            var_name = trim(dataEntry%short_name)//"_tavg"
            call PFIO_write_single_var(n, PFIOmodel_idx, vcol_id, TRIM(file_name), &
               dataEntry%modelOutput(1,:,1:nlev),  &
               var_name, local_var2D, local_var3D, local_var4D, &
               i1, i2, j1, j2, nlev, nmodel_status)

            CALL increment_field_counter()
            var_name2 = trim(dataEntry%short_name)//"_inst"
            call PFIO_write_single_var(n, PFIOmodel_idx, vcol_id, TRIM(file_name), &
               dataEntry%modelOutput(2,:,1:nlev),  &
               var_name2, local_var2D, local_var3D, local_var4D, &
               i1, i2, j1, j2, nlev, nmodel_status)

            ! time-averaged values or instantaneous values
         else
            CALL increment_field_counter(nmodel_status)
            call PFIO_write_single_var(n, PFIOmodel_idx, vcol_id, TRIM(file_name), &
               dataEntry%modelOutput(1,:,1:nlev),  &
               var_name, local_var2D, local_var3D, local_var4D, &
               i1, i2, j1, j2, nlev, nmodel_status)
         end if
         if (dataEntry%minmaxOpt.gt.0) then 
            CALL increment_field_counter()
            var_name2 = trim(dataEntry%short_name)//"_min"
            call PFIO_write_single_var(n, PFIOmodel_idx, vcol_id, TRIM(file_name), &
               dataEntry%minimum(:,1:nlev),  &
               var_name2, local_var2D, local_var3D, local_var4D, &
               i1, i2, j1, j2, nlev, nmodel_status)

            CALL increment_field_counter()
            var_name2 = trim(dataEntry%short_name)//"_max"
            call PFIO_write_single_var(n, PFIOmodel_idx, vcol_id, TRIM(file_name), &
               dataEntry%maximum(:,1:nlev),  &
               var_name2, local_var2D, local_var3D, local_var4D, &
               i1, i2, j1, j2, nlev, nmodel_status)
         endif
      end if IF_selectOpt

   end subroutine PFIO_write_variable

!EOC
! -----------------------------------------------------------------------
!BOP
! !ROUTINE: PFIO_write_single_var
! \label{PFIO_write_single_var}
! 
! !INTERFACE:
   subroutine PFIO_write_single_var(n, PFIOmodel_idx, vcol_id, file_name, var_data, var_name, &
         local_var2D, local_var3D, local_var4D, &
         i1, i2, j1, j2, nlev, nmodel_status)
! !USES: 

! !ARGUMENTS: 
      integer, intent(in) :: n
      integer, intent(in) :: PFIOmodel_idx
      integer, intent(in) :: vcol_id
      integer, intent(in) :: i1, i2, j1, j2
      integer, intent(in) :: nlev
      type(local_2Dfield_var), intent(inOut) :: local_var2D(:)
      type(local_3Dfield_var), intent(inOut) :: local_var3D(:)
      type(local_4Dfield_var), intent(inOut) :: local_var4D(:)
      character(len=*)    :: var_name
      character(len=*)    :: file_name
      real, intent(in)    :: var_data(LIS_rc%ntiles(n), nlev)
      integer, intent(in) :: nmodel_status
!
! !DESCRIPTION:
!  Write a real variable to a netcdf output file with some diagnostic 
!  statistics written to a text file. 
!
!  The arguments are: 
!  \begin{description}
!   \item [n]
!     index of the domain or nest.
!   \item [var_data]
!     variables being written, dimensioned in the tile space
!   \item[var_name]
!     name of the variable being written
!   \item [flag]
!    option to determine if the variable needs to be written (1-write, 
!    0-do not write)
!  \end{description}
!
! !LOCAL VARIABLES:
      integer              :: l, iret
      integer              :: global_dim(2)
      real                 :: vmean,vstdev,vmin,vmax
      real, allocatable    :: var1(:,:)
      real, allocatable    :: loclat(:)
      real, allocatable    :: loclon(:)
      real, allocatable    :: var1_ens(:,:,:)
      integer              :: count1 ,c,r,m,gid,ntiles,ierr,i,t
      integer              :: ews_ind, nss_ind
      integer              :: idx2, idx3, idx4, gindex
      type(ArrayReference) :: ref
!EOP
!------------------------------------------------------------------------------
!BOC

      ! Write output in 1-D ensemble tile array space:
      WOPT: if (LIS_rc%wopt.eq."1d tilespace") then !tiled output
         !iret = nf90_put_var(ftn,varid,gtmp1,(/1,dim1/), (/LIS_rc%glbntiles_red(n),1/))
         call map_1dtile_to_1darray(n, var_data, local_var2D(idx_field2d)%var2d(:,1:nlev), i1, i2, nlev)
         ref =  ArrayReference(local_var2D(idx_field2d)%var2d(:,1:nlev))
         call o_Clients%collective_stage_data(PFIO_bundle%hist_id(n,vcol_id,PFIOmodel_idx), &
            TRIM(file_name), TRIM(var_name), ref, &
            start        = [ i1, 1 ], &
            global_start = [  1, 1 ], & 
            global_count = [ LIS_rc%glbntiles_red(n), nlev ] )

         ! deallaocate(var1)
         ! Write output in 2d grid space:
      elseif(LIS_rc%wopt.eq."2d gridspace") then
         CALL get_LIS_globaldim(n, global_dim)
         allocate(var1(LIS_rc%ngrid(n),nlev))
         var1 = 0.0
         do i=1,LIS_rc%ntiles(n),LIS_rc%nensem(n)
            c = LIS_domain(n)%tile(i)%index
            do m=1,LIS_rc%nensem(n)
               t = i+m-1
               do l = 1, nlev
                  if ( var_data(t,l) == LIS_rc%udef) then
                     var1(c,l) = LIS_rc%udef ! pfio_missing_value 
                  else
                     var1(c,l) = var1(c,l) + var_data(t,l)*LIS_domain(n)%tile(t)%fgrd*LIS_domain(n)%tile(t)%pens
                  endif
               enddo
            enddo
         enddo

         ! The latlon fields are written to 1D
         if ((LIS_rc%nlatlon_dimensions == '1D') .AND. (nmodel_status > 0)) then
            if (nmodel_status.eq.1) then   ! lat
               allocate(loclat(LIS_rc%lnr(n)))
               loclat = LIS_rc%udef

               do r=1,LIS_rc%lnr(n)
                  do c=1,LIS_rc%lnc(n)
                     gindex = c+(r-1)*LIS_rc%lnc(n)
                     loclat(r) = LIS_domain(n)%lat(gindex)
                  enddo
               enddo

               ref = ArrayReference(loclat(:))
               call o_Clients%collective_stage_data(PFIO_bundle%hist_id(n,vcol_id,PFIOmodel_idx), &
                  TRIM(file_name), TRIM(var_name), ref, &
                  start        = [ j1 ], &
                  global_start = [ 1  ], & 
                  global_count = [ global_dim(2) ] )
               deallocate(loclat) 
            elseif(nmodel_status.eq.2) then !lon
               allocate(loclon(LIS_rc%lnc(n)))
               loclon = LIS_rc%udef

               do r=1,LIS_rc%lnr(n)
                  do c=1,LIS_rc%lnc(n)
                     gindex = c+(r-1)*LIS_rc%lnc(n)
                     loclon(c) = LIS_domain(n)%lon(gindex)
                  enddo
               enddo

               ref = ArrayReference(loclon(:))
               call o_Clients%collective_stage_data(PFIO_bundle%hist_id(n,vcol_id,PFIOmodel_idx), &
                  TRIM(file_name), TRIM(var_name), ref, &
                  start        = [ i1 ], &
                  global_start = [ 1  ], & 
                  global_count = [ global_dim(1) ] )
               deallocate(loclon) 
            endif
            ! The latlon fields are written to 2D
         else
            call map_1dtile_to_2darray(n, var1, local_var3D(idx_field3d)%var3d(:,:,1:nlev), &
               i1, i2, j1, j2, nlev)
            ref =  ArrayReference(local_var3D(idx_field3d)%var3d(:,:,1:nlev))

            call o_Clients%collective_stage_data(PFIO_bundle%hist_id(n,vcol_id,PFIOmodel_idx), &
               TRIM(file_name), TRIM(var_name), ref, &
               start        = [ i1, j1, 1 ], &
               global_start = [  1,  1, 1 ], & 
               global_count = [ global_dim(1), global_dim(2), nlev ] )

            deallocate(var1) 
         endif

         ! Write output in 2D ensemble grid space:
      elseif(LIS_rc%wopt.eq."2d ensemble gridspace") then
         CALL get_LIS_globaldim(n, global_dim)

         ! Non-model output field status (T=non-model; F=model-based):
         if (nmodel_status > 0) then   ! non-model output field status
            allocate(var1(LIS_rc%ngrid(n),nlev))
            var1 = 0.0
            do i=1,LIS_rc%ntiles(n),LIS_rc%nensem(n)
               c = LIS_domain(n)%tile(i)%index
               do m=1,LIS_rc%nensem(n)
                  t = i+m-1
                  do l = 1, nlev
                     if ( var_data(t,l) == LIS_rc%udef) then
                        var1(c,l) = LIS_rc%udef
                     else
                        var1(c,l) = var1(c,l) + var_data(t,l)*LIS_domain(n)%tile(t)%fgrd*LIS_domain(n)%tile(t)%pens
                     endif
                  enddo
               enddo
            enddo

            ! The latlon fields are written to 1D
            if (LIS_rc%nlatlon_dimensions == '1D') then
               if (nmodel_status.eq.1) then   ! lat
                  allocate(loclat(LIS_rc%lnr(n)))
                  loclat = LIS_rc%udef

                  do r=1,LIS_rc%lnr(n)
                     do c=1,LIS_rc%lnc(n)
                        gindex = c+(r-1)*LIS_rc%lnc(n)
                        loclat(r) = LIS_domain(n)%lat(gindex)
                     enddo
                  enddo

                  ref = ArrayReference(loclat(:))
                  call o_Clients%collective_stage_data(PFIO_bundle%hist_id(n,vcol_id,PFIOmodel_idx), &
                     TRIM(file_name), TRIM(var_name), ref, &
                     start        = [ j1 ], &
                     global_start = [ 1  ], &
                     global_count = [ global_dim(2) ] )
                  deallocate(loclat)
               elseif(nmodel_status.eq.2) then !lon
                  allocate(loclon(LIS_rc%lnc(n)))
                  loclon = LIS_rc%udef

                  do r=1,LIS_rc%lnr(n)
                     do c=1,LIS_rc%lnc(n)
                        gindex = c+(r-1)*LIS_rc%lnc(n)
                        loclon(c) = LIS_domain(n)%lon(gindex)
                     enddo
                  enddo

                  ref = ArrayReference(loclon(:))
                  call o_Clients%collective_stage_data(PFIO_bundle%hist_id(n,vcol_id,PFIOmodel_idx), &
                     TRIM(file_name), TRIM(var_name), ref, &
                     start        = [ i1 ], &
                     global_start = [ 1  ], &
                     global_count = [ global_dim(1) ] )
                  deallocate(loclon)
               endif
               ! The latlon fields are written to 2D
            else
               call map_1dtile_to_2darray(n, var1, local_var3D(idx_field3d)%var3d(:,:,1:nlev), i1, i2, j1, j2, nlev)
               ref =  ArrayReference(local_var3D(idx_field3d)%var3d(:,:,1:nlev))

               call o_Clients%collective_stage_data(PFIO_bundle%hist_id(n,vcol_id,PFIOmodel_idx), &
                  TRIM(file_name), TRIM(var_name), ref, &
                  start        = [ i1, j1, 1 ], &
                  global_start = [  1,  1, 1 ], &
                  global_count = [ global_dim(1), global_dim(2), nlev ] )

               deallocate(var1)
            endif
            ! Model-based field output:
         else
            allocate(var1_ens(LIS_rc%ngrid(n), LIS_rc%nensem(n), nlev))

            var1_ens = 0
            do i=1,LIS_rc%ntiles(n),LIS_rc%nensem(n)
               c = LIS_domain(n)%tile(i)%index
               do m=1,LIS_rc%nensem(n)
                  t = i+m-1
                  do l =1, nlev
                     if ( var_data(t,l) == LIS_rc%udef ) then
                        var1_ens(c,m,l) = LIS_rc%udef
                     else
                        var1_ens(c,m,l) =  var_data(t,l)*LIS_domain(n)%tile(t)%fgrd
                     endif
                  enddo
               enddo
            enddo

            do m=1,LIS_rc%nensem(n)
               call map_1dtile_to_2darray(n, var1_ens(:,m,:), local_var4D(idx_field4d)%var4d(:,:,m,1:nlev), i1, i2, j1, j2, nlev)
            enddo

            m = LIS_rc%nensem(n)
            ref =  ArrayReference(local_var4D(idx_field4d)%var4d(:,:,:,1:nlev))
            call o_Clients%collective_stage_data(PFIO_bundle%hist_id(n,vcol_id,PFIOmodel_idx), &
               TRIM(file_name), TRIM(var_name), ref, &
               start        = [ i1, j1, 1, 1 ], &
               global_start = [  1,  1, 1, 1 ], & 
               global_count = [ global_dim(1), global_dim(2), m, nlev] )

            deallocate(var1_ens) 
         endif

      end if WOPT


   end subroutine PFIO_write_single_var
!EOC
!------------------------------------------------------------------------------
!BOP
! !ROUTINE: LIS_PFIO_write_routingdata
! \label{LIS_PFIO_write_routingdata}
!
! !INTERFACE: 
   subroutine LIS_PFIO_write_routingdata(n, PFIOmodel_idx, vcol_id, file_name, outInterval, group)

      integer,          intent(in) :: n, PFIOmodel_idx
      integer,          intent(in) :: vcol_id 
      character(len=*), intent(in) :: file_name
      real,             intent(in) :: outInterval
      integer,          intent(in) :: group
! 
! !DESCRIPTION: 
!  This routine writes variables to a NETCDF file
!
! !LOCAL VARIABLES:
      real                    :: time_data(1)
      integer                 :: i,k,m,t, group_, status, rc
      character(len=8)        :: xtime_begin_date
      character(len=6)        :: xtime_begin_time
      character(len=50)       :: xtime_units
      character(len=50)       :: xtime_timeInc
      character(len=8)        :: date
      character(len=10)       :: time
      character(len=5)        :: zone
      integer, dimension(8)   :: values
      Type(Variable)          :: v
      type(StringVariableMap) :: var_map
      type(LIS_metadataEntry), pointer :: dataEntry
      type(local_2Dfield_var), allocatable :: local_var2D(:)
      type(local_3Dfield_var), allocatable :: local_var3D(:)
      type(local_4Dfield_var), allocatable :: local_var4D(:)
      integer                 :: i1, i2, j1, j2
      type(ArrayReference)    :: ref
      integer                 :: global_dim(2)
      integer                 :: r, c, index1
!EOP
!---------------------------------------------------------------------------------------------
!BOC
      write(LIS_logunit,'(a,i,a,a)')'[INFO-PFIO] Writing routing output [',PFIOmodel_idx,'] to:  ', TRIM(file_name)

      ! Update the time variable
      !-------------------------
      call date_and_time(date, time, zone, values)
      write(xtime_units, 200) LIS_rc%yr, LIS_rc%mo, LIS_rc%da, LIS_rc%hr, LIS_rc%mn, LIS_rc%ss
      200   format ('minutes since ',I4.4,'-',I2.2,'-',I2.2,' ',I2.2,':', I2.2,':',I2.2)
      write(xtime_begin_date, fmt='(I4.4,I2.2,I2.2)') LIS_rc%yr, LIS_rc%mo, LIS_rc%da
      write(xtime_begin_time, fmt='(I2.2,I2.2,I2.2)') LIS_rc%hr, LIS_rc%mn, LIS_rc%ss
      write(xtime_timeInc, fmt='(I20)')  nint(outInterval)

      time_data = (/ 1.0 /)
      v = Variable(type=PFIO_REAL32, dimensions='time')
      call v%add_attribute("units",trim(xtime_units))
      call v%add_attribute("long_name","time")
      call v%add_attribute("time_increment",trim(adjustl(xtime_timeInc)))
      call v%add_attribute("begin_date",xtime_begin_date)
      call v%add_attribute("begin_time",xtime_begin_time)
      call v%add_const_value(UnlimitedEntity( time_data ))
      call PFIO_bundle%fmd(n,vcol_id,PFIOmodel_idx)%add_variable('time', v)
      call var_map%insert('time', v)
      call o_Clients%modify_metadata(PFIO_bundle%hist_id(n,vcol_id,PFIOmodel_idx), var_map=var_map, rc=status)

      call o_Clients%set_optimal_server(nwriting=1)

      group_ = group

      if (LIS_rc%wopt.eq."1d tilespace") then
         CALL get_LIS_interior_tile(n, i1, i2)
      elseif ( (LIS_rc%wopt.eq."2d gridspace") .OR. (LIS_rc%wopt.eq."2d ensemble gridspace") )then 
         CALL get_LIS_interior_grid(n, i1, i2, j1, j2)
      endif
      CALL get_LIS_globaldim(n, global_dim)


      dataEntry => LIS_histData(n)%head_routing_list

      SELECT CASE(LIS_rc%wopt)
      CASE("2d gridspace")
         ALLOCATE(local_var3D(total_num_3Dfields(PFIOmodel_idx)))
      CASE("2d ensemble gridspace")
         ALLOCATE(local_var3D(total_num_3Dfields(PFIOmodel_idx)))
         ALLOCATE(local_var4D(total_num_4Dfields(PFIOmodel_idx)))
      END SELECT
      call allocate_local_var(n, PFIOmodel_idx, local_var2D, local_var3D, local_var4D, i1, i2, j1, j2)
      idx_field2d = 0
      idx_field3d = 0
      idx_field4d = 0

      ! Writing latitude/longitude
      !---------------------------
      !if ( (LIS_rc%wopt.eq."2d gridspace") .OR. (LIS_rc%wopt.eq."2d ensemble gridspace") )then 
      call PFIO_write_variable(n, PFIOmodel_idx, vcol_id, xlat,  file_name, local_var2D, local_var3D, local_var4D, &
         i1, i2, j1, j2, non_model_fields=1)
      call PFIO_write_variable(n, PFIOmodel_idx, vcol_id, xlong, file_name, local_var2D, local_var3D, local_var4D, &
         i1, i2, j1, j2, non_model_fields=2)
      !endif

      ! Writing the model fields
      !-------------------------
      do while ( associated(dataEntry) )
         call PFIO_write_routingvariable(n, PFIOmodel_idx, vcol_id, dataEntry, file_name, &
            local_var2D, local_var3D, local_var4D, i1, i2, j1, j2)
         dataEntry => dataEntry%next
      enddo

      ! write in the file and close it
      call o_Clients%done_collective_stage()

      call o_Clients%post_wait()

      ! After writing reset the variables
      SELECT CASE(LIS_rc%wopt)
      CASE("1d tilespace")
         DO i= 1, total_num_2Dfields(PFIOmodel_idx)
            deallocate(local_var2D(i)%var2d)
         ENDDO
         deallocate(local_var2D)
      CASE("2d gridspace")
         DO i= 1, total_num_3Dfields(PFIOmodel_idx)
            deallocate(local_var3D(i)%var3d)
         ENDDO
         deallocate(local_var3D)
      CASE("2d ensemble gridspace")
         DO i= 1, total_num_3Dfields(PFIOmodel_idx)
            deallocate(local_var3D(i)%var3d)
         ENDDO
         deallocate(local_var3D)
         DO i= 1, total_num_4Dfields(PFIOmodel_idx)
            deallocate(local_var4D(i)%var4d)
         ENDDO
         deallocate(local_var4D)
      END SELECT

      call LIS_resetOutputVars(n, group_)

   end subroutine LIS_PFIO_write_routingdata
!EOC
!------------------------------------------------------------------------------
!BOP
! !ROUTINE: PFIO_write_routingvariable
! \label{PFIO_write_routingvariable}
!
! !INTERFACE: 
   subroutine PFIO_write_routingvariable(n, PFIOmodel_idx, vcol_id, dataEntry, file_name, &
         local_var2D, local_var3D, local_var4D, &
         i1, i2, j1, j2, non_model_fields)

      integer,   intent(in)   :: n 
      integer,   intent(in)   :: PFIOmodel_idx 
      integer,   intent(in)   :: vcol_id 
      integer,   intent(in)   :: i1, i2, j1, j2
      type(local_2Dfield_var), intent(inOut) :: local_var2D(:)
      type(local_3Dfield_var), intent(inOut) :: local_var3D(:)
      type(local_4Dfield_var), intent(inOut) :: local_var4D(:)
      type(LIS_metadataEntry), pointer :: dataEntry
      character(len=*), intent(in) :: file_name
      integer, optional, intent(in) :: non_model_fields
! 
! !DESCRIPTION: 
!  This routine writes a single variable to a NETCDF file
!  The arguments are: 
!  \begin{description}
!    \item[n] index of the nest
!   \item[dataEntry]
!    object containing the values and attributes of the variable to be 
!    written
!  \end{description}
!
!   The routines invoked are: 
!   \begin{description}
!   \item[LIS\_writevar\_netcdf](\ref{LIS_writevar_netcdf})
!     writes a variable into a netcdf formatted file. 
!   \end{description}
!
! !LOCAL VARIABLES:
      integer       :: i,k,m,t, nlev
      character(len=90) :: var_name
      character(len=90) :: var_name2
      integer       :: nmodel_status
!EOP    
!------------------------------------------------------------------------------
!BOC

      nmodel_status = 0
      if (present(non_model_fields)) nmodel_status = non_model_fields

      IF_selectOpt: if (dataEntry%selectOpt.eq.1) then

         IF (nmodel_status == 0) THEN
            do t=1,LIS_rc%nroutinggrid(n)*LIS_rc%nensem(n)
               do k=1,dataEntry%vlevels
                  if (dataEntry%count(t,k).gt.0) then 
                     SELECT CASE(dataEntry%timeAvgOpt)
                     CASE(3)
                        continue   
                     CASE (1, 2)
                        dataEntry%modelOutput(1,t,k) = dataEntry%modelOutput(1,t,k)/ dataEntry%count(t,k)
                     CASE DEFAULT
                        continue   
                     END SELECT
                  else
                     dataEntry%modelOutput(1,t,k) = pfio_missing_value
                  endif
               enddo
            enddo
         ENDIF

         if (nmodel_status == 0) then
            SELECT CASE(dataEntry%timeAvgOpt)
            CASE(0)
               var_name = trim(dataEntry%short_name)//'_inst'
            CASE(1)
               var_name = trim(dataEntry%short_name)//'_tavg'
            CASE(3)
               var_name = trim(dataEntry%short_name)//'_acc'
            END SELECT
         else
            var_name = trim(dataEntry%short_name)
         endif

         ! accumulated values
         ! time-averaged values and instantaneous values
         nlev = dataEntry%vlevels
         if (dataEntry%timeAvgOpt.eq.2) then 
            CALL increment_field_counter()
            var_name = trim(dataEntry%short_name)//"_tavg"
            call PFIO_write_single_routingvar(n, PFIOmodel_idx, vcol_id, TRIM(file_name), &
               dataEntry%modelOutput(1,:,1:nlev),  &
               var_name, local_var2D, local_var3D, local_var4D, &
               i1, i2, j1, j2, nlev, nmodel_status)

            CALL increment_field_counter()
            var_name2 = trim(dataEntry%short_name)//"_inst"
            call PFIO_write_single_routingvar(n, PFIOmodel_idx, vcol_id, TRIM(file_name), &
               dataEntry%modelOutput(2,:,1:nlev),  &
               var_name2, local_var2D, local_var3D, local_var4D, &
               i1, i2, j1, j2, nlev, nmodel_status)

            ! time-averaged values or instantaneous values
         else
            CALL increment_field_counter(nmodel_status)
            call PFIO_write_single_routingvar(n, PFIOmodel_idx, vcol_id, TRIM(file_name), &
               dataEntry%modelOutput(1,:,1:nlev),  &
               var_name, local_var2D, local_var3D, local_var4D, &
               i1, i2, j1, j2, nlev, nmodel_status)
         end if
         if (dataEntry%minmaxOpt.gt.0) then 
            CALL increment_field_counter()
            var_name2 = trim(dataEntry%short_name)//"_min"
            call PFIO_write_single_routingvar(n, PFIOmodel_idx, vcol_id, TRIM(file_name), &
               dataEntry%minimum(:,1:nlev),  &
               var_name2, local_var2D, local_var3D, local_var4D, &
               i1, i2, j1, j2, nlev, nmodel_status)

            CALL increment_field_counter()
            var_name2 = trim(dataEntry%short_name)//"_max"
            call PFIO_write_single_routingvar(n, PFIOmodel_idx, vcol_id, TRIM(file_name), &
               dataEntry%maximum(:,1:nlev),  &
               var_name2, local_var2D, local_var3D, local_var4D, &
               i1, i2, j1, j2, nlev, nmodel_status)
         endif
      end if IF_selectOpt

   end subroutine PFIO_write_routingvariable
!EOC
! -----------------------------------------------------------------------
!BOP
! !ROUTINE: PFIO_write_single_var
! \label{PFIO_write_single_var}
! 
! !INTERFACE:
   subroutine PFIO_write_single_routingvar(n, PFIOmodel_idx, vcol_id, file_name, var_data, var_name, &
         local_var2D, local_var3D, local_var4D, &
         i1, i2, j1, j2, nlev, nmodel_status)
! !USES: 

! !ARGUMENTS: 
      integer, intent(in) :: n, PFIOmodel_idx
      integer, intent(in) :: vcol_id
      integer, intent(in) :: i1, i2, j1, j2
      integer, intent(in) :: nlev
      type(local_2Dfield_var), intent(inOut) :: local_var2D(:)
      type(local_3Dfield_var), intent(inOut) :: local_var3D(:)
      type(local_4Dfield_var), intent(inOut) :: local_var4D(:)
      character(len=*)    :: var_name
      character(len=*)    :: file_name
      real, intent(in)    :: var_data(LIS_rc%ntiles(n), nlev)
      integer, intent(in) :: nmodel_status
!
! !DESCRIPTION:
!  Write a real variable to a netcdf output file with some diagnostic 
!  statistics written to a text file. 
!
!  The arguments are: 
!  \begin{description}
!   \item [n]
!     index of the domain or nest.
!   \item [var_data]
!     variables being written, dimensioned in the tile space
!   \item[var_name]
!     name of the variable being written
!   \item [flag]
!    option to determine if the variable needs to be written (1-write, 
!    0-do not write)
!  \end{description}
!
! !LOCAL VARIABLES:
      integer              :: l, iret
      integer              :: global_dim(2)
      real                 :: vmean,vstdev,vmin,vmax
      real, allocatable    :: var1(:,:)
      real, allocatable    :: loclat(:)
      real, allocatable    :: loclon(:)
      real, allocatable    :: var1_ens(:,:,:)
      integer              :: count1 ,c,r,m,gid,ntiles,ierr,i,t
      integer              :: ews_ind, nss_ind
      integer              :: idx2, idx3, idx4, gindex
      type(ArrayReference) :: ref
!EOP
!------------------------------------------------------------------------------
!BOC

      ! Write output in 2d grid space:
      WOPT: if (LIS_rc%wopt.eq."2d gridspace") then
         CALL get_LIS_globaldim(n, global_dim)
         allocate(var1(LIS_rc%nroutinggrid(n),nlev))
         var1 = 0.0
         do i=1,LIS_rc%nroutinggrid(n)
            do m=1,LIS_rc%nensem(n)
               t = m+(i-1)*LIS_rc%nensem(n)
               do l = 1, nlev
                  if ( var_data(t,l) == LIS_rc%udef) then
                     var1(i,l) = LIS_rc%udef ! pfio_missing_value 
                  else
                     var1(i,l) = var1(i,l) + var_data(t,l)*LIS_routing(n)%tile(t)%fgrd*LIS_routing(n)%tile(t)%pens
                  endif
               enddo
            enddo
         enddo

         ! The latlon fields are written to 1D
         if ((LIS_rc%nlatlon_dimensions == '1D') .AND. (nmodel_status > 0)) then
            if (nmodel_status.eq.1) then   ! lat
               allocate(loclat(LIS_rc%lnr(n)))
               loclat = LIS_rc%udef

               do r=1,LIS_rc%lnr(n)
                  do c=1,LIS_rc%lnc(n)
                     gindex = c+(r-1)*LIS_rc%lnc(n)
                     loclat(r) = LIS_domain(n)%lat(gindex)
                  enddo
               enddo

               ref = ArrayReference(loclat(:))
               call o_Clients%collective_stage_data(PFIO_bundle%hist_id(n,vcol_id,PFIOmodel_idx), &
                  TRIM(file_name), TRIM(var_name), ref, &
                  start        = [ j1 ], &
                  global_start = [ 1  ], & 
                  global_count = [ global_dim(2) ] )
               deallocate(loclat) 
            elseif(nmodel_status.eq.2) then !lon
               allocate(loclon(LIS_rc%lnc(n)))
               loclon = LIS_rc%udef

               do r=1,LIS_rc%lnr(n)
                  do c=1,LIS_rc%lnc(n)
                     gindex = c+(r-1)*LIS_rc%lnc(n)
                     loclon(c) = LIS_domain(n)%lon(gindex)
                  enddo
               enddo

               ref = ArrayReference(loclon(:))
               call o_Clients%collective_stage_data(PFIO_bundle%hist_id(n,vcol_id,PFIOmodel_idx), &
                  TRIM(file_name), TRIM(var_name), ref, &
                  start        = [ i1 ], &
                  global_start = [ 1  ], & 
                  global_count = [ global_dim(1) ] )
               deallocate(loclon) 
            endif
            ! The latlon fields are written to 2D
         else
            call map_1dtile_to_2darray(n, var1, local_var3D(idx_field3d)%var3d(:,:,1:nlev), &
               i1, i2, j1, j2, nlev)
            ref =  ArrayReference(local_var3D(idx_field3d)%var3d(:,:,1:nlev))

            call o_Clients%collective_stage_data(PFIO_bundle%hist_id(n,vcol_id,PFIOmodel_idx), &
               TRIM(file_name), TRIM(var_name), ref, &
               start        = [ i1, j1, 1 ], &
               global_start = [  1,  1, 1 ], & 
               global_count = [ global_dim(1), global_dim(2), nlev ] )

            deallocate(var1) 
         endif

         ! Write output in 2D ensemble grid space:
      elseif(LIS_rc%wopt.eq."2d ensemble gridspace") then
         CALL get_LIS_globaldim(n, global_dim)

         ! Non-model output field status (T=non-model; F=model-based):
         if (nmodel_status > 0) then   ! non-model output field status
            allocate(var1(LIS_rc%nroutinggrid(n),nlev))
            var1 = 0.0
            do i=1,LIS_rc%nroutinggrid(n)
               do m=1,LIS_rc%nensem(n)
                  t = m + (i-1)*LIS_rc%nensem(n)
                  do l = 1, nlev
                     if ( var_data(t,l) == LIS_rc%udef) then
                        var1(i,l) = LIS_rc%udef
                     else
                        var1(i,l) = var1(i,l) + var_data(t,l)*LIS_routing(n)%tile(t)%fgrd*LIS_routing(n)%tile(t)%pens
                     endif
                  enddo
               enddo
            enddo

            ! The latlon fields are written to 1D
            if (LIS_rc%nlatlon_dimensions == '1D') then
               if (nmodel_status.eq.1) then   ! lat
                  allocate(loclat(LIS_rc%lnr(n)))
                  loclat = LIS_rc%udef

                  do r=1,LIS_rc%lnr(n)
                     do c=1,LIS_rc%lnc(n)
                        gindex = c+(r-1)*LIS_rc%lnc(n)
                        loclat(r) = LIS_domain(n)%lat(gindex)
                     enddo
                  enddo

                  ref = ArrayReference(loclat(:))
                  call o_Clients%collective_stage_data(PFIO_bundle%hist_id(n,vcol_id,PFIOmodel_idx), &
                     TRIM(file_name), TRIM(var_name), ref, &
                     start        = [ j1 ], &
                     global_start = [ 1  ], &
                     global_count = [ global_dim(2) ] )
                  deallocate(loclat)
               elseif(nmodel_status.eq.2) then !lon
                  allocate(loclon(LIS_rc%lnc(n)))
                  loclon = LIS_rc%udef

                  do r=1,LIS_rc%lnr(n)
                     do c=1,LIS_rc%lnc(n)
                        gindex = c+(r-1)*LIS_rc%lnc(n)
                        loclon(c) = LIS_domain(n)%lon(gindex)
                     enddo
                  enddo

                  ref = ArrayReference(loclon(:))
                  call o_Clients%collective_stage_data(PFIO_bundle%hist_id(n,vcol_id,PFIOmodel_idx), &
                     TRIM(file_name), TRIM(var_name), ref, &
                     start        = [ i1 ], &
                     global_start = [ 1  ], &
                     global_count = [ global_dim(1) ] )
                  deallocate(loclon)
               endif
               ! The latlon fields are written to 2D
            else
               call map_1dtile_to_2darray(n, var1, local_var3D(idx_field3d)%var3d(:,:,1:nlev), i1, i2, j1, j2, nlev)
               ref =  ArrayReference(local_var3D(idx_field3d)%var3d(:,:,1:nlev))

               call o_Clients%collective_stage_data(PFIO_bundle%hist_id(n,vcol_id,PFIOmodel_idx), &
                  TRIM(file_name), TRIM(var_name), ref, &
                  start        = [ i1, j1, 1 ], &
                  global_start = [  1,  1, 1 ], &
                  global_count = [ global_dim(1), global_dim(2), nlev ] )

               deallocate(var1)
            endif
            ! Model-based field output:
         else
            allocate(var1_ens(LIS_rc%ngrid(n), LIS_rc%nensem(n), nlev))

            var1_ens = 0
            do i=1,LIS_rc%ntiles(n),LIS_rc%nensem(n)
               c = LIS_domain(n)%tile(i)%index
               do m=1,LIS_rc%nensem(n)
                  t = i+m-1
                  do l =1, nlev
                     if ( var_data(t,l) == LIS_rc%udef ) then
                        var1_ens(c,m,l) = LIS_rc%udef
                     else
                        var1_ens(c,m,l) =  var_data(t,l)*LIS_domain(n)%tile(t)%fgrd
                     endif
                  enddo
               enddo
            enddo

            do m=1,LIS_rc%nensem(n)
               call map_1dtile_to_2darray(n, var1_ens(:,m,:), local_var4D(idx_field4d)%var4d(:,:,m,1:nlev), i1, i2, j1, j2, nlev)
            enddo

            m = LIS_rc%nensem(n)
            ref =  ArrayReference(local_var4D(idx_field4d)%var4d(:,:,:,1:nlev))
            call o_Clients%collective_stage_data(PFIO_bundle%hist_id(n,vcol_id,PFIOmodel_idx), &
               TRIM(file_name), TRIM(var_name), ref, &
               start        = [ i1, j1, 1, 1 ], &
               global_start = [  1,  1, 1, 1 ], & 
               global_count = [ global_dim(1), global_dim(2), m, nlev] )

            deallocate(var1_ens) 
         endif

      end if WOPT

   end subroutine PFIO_write_single_routingvar
!EOC
! -----------------------------------------------------------------------
#if 0
! LIS' restart format is not compatible with PFIO.
! Save these routines for future reference.
!
!BOP
! !ROUTINE: PFIO_restart_metadata
! \label{PFIO_restart_metadata}
! 
! !INTERFACE: PFIO_create_restart_metadata
   subroutine PFIO_restart_metadata(n, m, model_name, var_flag, &
         dim1, dim2, dim3, dim4, dim5, dim6, dim7, dim8, dim9, dim10,&
         output_format)
! !USES: 
!
! !INPUTS PARAMETERS: 
      integer,   intent(in)     :: n
      integer,   intent(in)     :: m
      character(len=*), intent(in) :: model_name
      character(len=*), optional :: var_flag
      integer,    optional      :: dim1
      integer,    optional      :: dim2
      integer,    optional      :: dim3
      integer,    optional      :: dim4
      integer,    optional      :: dim5
      integer,    optional      :: dim6
      integer,    optional      :: dim7
      integer,    optional      :: dim8
      integer,    optional      :: dim9
      integer,    optional      :: dim10
      character(len=*), optional :: output_format

! 
! !DESCRIPTION: 
!  This routine creates a PFIO object that contains metadata to be used
!  in a netCDF restart file.
!  \begin{description}
!    \item[n] index of the nest
!    \item[ftn] file unit for the output file
!    \item[model\_name] Name of the model that generates the output
!  \end{description}
!
!   The routines invoked are: 
!   \begin{description}
!   \item[defineNETCDFheadervar](\ref{defineNETCDFheaderVar})
!     writes the required headers for a single variable
!   \item[writeSingleNETCDFvar](\ref{writeSingleNETCDFvar})
!     writes a single variable into a netcdf formatted file. 
!   \item[LIS\_verify](\ref{LIS_verify})
!     call to check if the return value is valid or not.
!   \end{description}
!
! !LOCAL VARIABLES:
      Type(Variable)          :: v

      integer, dimension(8) :: values
      character(len=8)  :: date
      character(len=10) :: time
      character(len=5)  :: zone
      character*20      :: wout
      character*50 :: var_flag_tmp

      real                    :: time_data(1)
      REAL                    :: SOUTH_WEST_CORNER_LAT
      REAL                    :: SOUTH_WEST_CORNER_LON
      integer                 :: ierr, idx, status
      integer                 :: deflate_level
      character(len=8)        :: xtime_begin_date
      character(len=6)        :: xtime_begin_time
      character(len=50)       :: xtime_units
      character(len=50)       :: xtime_timeInc
!EOP
!---------------------------------------------------------------------------------------------
!BOC

      if (present(var_flag)) then
         var_flag_tmp = var_flag
      else
         var_flag_tmp = ""
      endif


      if(.NOT.PRESENT(output_format)) then
         wout = LIS_rc%wout
      else
         wout = trim(output_format)
      endif

      call PFIO_bundle%fmd_rst(n)%add_dimension('ntiles', LIS_rc%glbnpatch_red(n,m), rc=status)

      if (present(dim1))  call PFIO_bundle%fmd_rst(n)%add_dimension('dim1',  dim1,  rc=status)
      if (present(dim2))  call PFIO_bundle%fmd_rst(n)%add_dimension('dim2',  dim2,  rc=status)
      if (present(dim3))  call PFIO_bundle%fmd_rst(n)%add_dimension('dim3',  dim3,  rc=status)
      if (present(dim4))  call PFIO_bundle%fmd_rst(n)%add_dimension('dim4',  dim4,  rc=status)
      if (present(dim5))  call PFIO_bundle%fmd_rst(n)%add_dimension('dim5',  dim5,  rc=status)
      if (present(dim6))  call PFIO_bundle%fmd_rst(n)%add_dimension('dim6',  dim6,  rc=status)
      if (present(dim7))  call PFIO_bundle%fmd_rst(n)%add_dimension('dim7',  dim7,  rc=status)
      if (present(dim8))  call PFIO_bundle%fmd_rst(n)%add_dimension('dim8',  dim8,  rc=status)
      if (present(dim9))  call PFIO_bundle%fmd_rst(n)%add_dimension('dim9',  dim9,  rc=status)
      if (present(dim10)) call PFIO_bundle%fmd_rst(n)%add_dimension('dim10', dim10, rc=status)

      deflate_level = LIS_rc%nc_deflate_lvl

      v = Variable(type=PFIO_REAL32, dimensions='ntiles', deflation=deflate_level)
      call v%add_attribute('units', 'degree_north')
      call v%add_attribute('long_name', 'latitude')
      call v%add_attribute('standard_name', 'latitude')
      call v%add_attribute('scale_factor', 1.0)
      call v%add_attribute('add_offset', 0.0)
      call v%add_attribute('missing_value', LIS_rc%udef)
      call v%add_attribute('_FillValue', LIS_rc%udef)
      call v%add_attribute('vmin', pfio_vmin)
      call v%add_attribute('vmax', pfio_vmax)
      call PFIO_bundle%fmd_rst(n)%add_variable('lat', v)

      v = Variable(type=PFIO_REAL32, dimensions='ntiles', deflation=deflate_level)
      call v%add_attribute('units', 'degree_east')
      call v%add_attribute('long_name', 'longitude')
      call v%add_attribute('standard_name', 'longitude')
      call v%add_attribute('scale_factor', 1.0)
      call v%add_attribute('add_offset', 0.0)
      call v%add_attribute('missing_value', LIS_rc%udef)
      call v%add_attribute('_FillValue', LIS_rc%udef)
      call v%add_attribute('vmin', pfio_vmin)
      call v%add_attribute('vmax', pfio_vmax)
      call PFIO_bundle%fmd_rst(n)%add_variable('lon', v)

      ! defining time field
      call date_and_time(date,time,zone,values)
      write(xtime_units, 200) LIS_rc%yr, LIS_rc%mo, LIS_rc%da, LIS_rc%hr, LIS_rc%mn, LIS_rc%ss
      200   format ('minutes since ',I4.4,'-',I2.2,'-',I2.2,' ',I2.2,':', I2.2,':',I2.2)
      write(xtime_begin_date, fmt='(I4.4,I2.2,I2.2)') LIS_rc%yr, LIS_rc%mo, LIS_rc%da
      write(xtime_begin_time, fmt='(I2.2,I2.2,I2.2)') LIS_rc%hr, LIS_rc%mn, LIS_rc%ss
      !write(xtime_timeInc, fmt='(I20)')  nint(outInterval)

      time_data = 0.0
      v = Variable(type=PFIO_REAL32, dimensions='time')
      call v%add_attribute("units", trim(xtime_units))
      call v%add_attribute("long_name", "time")
      !call v%add_attribute("time_increment", trim(adjustl(xtime_timeInc)))
      call v%add_attribute("begin_date", xtime_begin_date)
      call v%add_attribute("begin_time", xtime_begin_time)
      call v%add_const_value(UnlimitedEntity( time_data ))
      call PFIO_bundle%fmd_rst(n)%add_variable('time', v)

      !------------------
      ! Global attributes
      !------------------
      ! The history attribute might be updated when the actual file is created.
      call PFIO_bundle%fmd_rst(n)%add_attribute("missing_value", LIS_rc%udef)
      call PFIO_bundle%fmd_rst(n)%add_attribute("title", "LIS land surface model restart")
      call PFIO_bundle%fmd_rst(n)%add_attribute("institution", trim(LIS_rc%institution))
      call PFIO_bundle%fmd_rst(n)%add_attribute("source", trim(model_name))
      call PFIO_bundle%fmd_rst(n)%add_attribute("history",  "created on date: "//date(1:4)//"-"//date(5:6)//"-"// date(7:8)//"T"//time(1:2)//":"//time(3:4)//":"//time(5:10))
      call PFIO_bundle%fmd_rst(n)%add_attribute("references", "Kumar_etal_EMS_2006, Peters-Lidard_etal_ISSE_2007")
      call PFIO_bundle%fmd_rst(n)%add_attribute("conventions", "CF-1.6")
      call PFIO_bundle%fmd_rst(n)%add_attribute("comment", "website: http://lis.gsfc.nasa.gov/")

      !============================================================================
      ! Make sure that all the processes have the same value the lower left corner.
      ! Broadcast the the root process values to all the other processes.
      ! This is only needed to write global attributes in the file.
      !============================================================================
      IF (LIS_masterproc) THEN
         SOUTH_WEST_CORNER_LAT = LIS_rc%gridDesc(n,4)
         SOUTH_WEST_CORNER_LON = LIS_rc%gridDesc(n,5)
      ELSE
         SOUTH_WEST_CORNER_LAT = -9999.0
         SOUTH_WEST_CORNER_LON = -9999.0
      ENDIF
      call MPI_Bcast(SOUTH_WEST_CORNER_LAT,1, MPI_REAL, 0, LIS_mpi_comm, ierr)
      call MPI_Bcast(SOUTH_WEST_CORNER_LON,1, MPI_REAL, 0, LIS_mpi_comm, ierr)

      if (LIS_rc%lis_map_proj.eq."latlon") then   ! latlon
         call PFIO_bundle%fmd_rst(n)%add_attribute("MAP_PROJECTION", "EQUIDISTANT CYLINDRICAL")
         call PFIO_bundle%fmd_rst(n)%add_attribute("SOUTH_WEST_CORNER_LAT", SOUTH_WEST_CORNER_LAT)
         call PFIO_bundle%fmd_rst(n)%add_attribute("SOUTH_WEST_CORNER_LON", SOUTH_WEST_CORNER_LON)
         call PFIO_bundle%fmd_rst(n)%add_attribute("DX", LIS_rc%gridDesc(n,9))
         call PFIO_bundle%fmd_rst(n)%add_attribute("DY", LIS_rc%gridDesc(n,10))       
      elseif (LIS_rc%lis_map_proj.eq."mercator") then 
         call PFIO_bundle%fmd_rst(n)%add_attribute("MAP_PROJECTION", "MERCATOR")
         call PFIO_bundle%fmd_rst(n)%add_attribute("SOUTH_WEST_CORNER_LAT", SOUTH_WEST_CORNER_LAT)
         call PFIO_bundle%fmd_rst(n)%add_attribute("SOUTH_WEST_CORNER_LON", SOUTH_WEST_CORNER_LON)
         call PFIO_bundle%fmd_rst(n)%add_attribute("TRUELAT1", LIS_rc%gridDesc(n,10))
         call PFIO_bundle%fmd_rst(n)%add_attribute("STANDARD_LON", LIS_rc%gridDesc(n,11))
         call PFIO_bundle%fmd_rst(n)%add_attribute("DX", LIS_rc%gridDesc(n,8))
         call PFIO_bundle%fmd_rst(n)%add_attribute("DY", LIS_rc%gridDesc(n,9))
      elseif (LIS_rc%lis_map_proj.eq."lambert") then ! lambert conformal
         call PFIO_bundle%fmd_rst(n)%add_attribute("MAP_PROJECTION",  "LAMBERT CONFORMAL")
         call PFIO_bundle%fmd_rst(n)%add_attribute("SOUTH_WEST_CORNER_LAT", SOUTH_WEST_CORNER_LAT)
         call PFIO_bundle%fmd_rst(n)%add_attribute("SOUTH_WEST_CORNER_LON", SOUTH_WEST_CORNER_LON)
         call PFIO_bundle%fmd_rst(n)%add_attribute("TRUELAT1", LIS_rc%gridDesc(n,10))
         call PFIO_bundle%fmd_rst(n)%add_attribute("TRUELAT2", LIS_rc%gridDesc(n,7))
         call PFIO_bundle%fmd_rst(n)%add_attribute("STANDARD_LON", LIS_rc%gridDesc(n,11))
         call PFIO_bundle%fmd_rst(n)%add_attribute("DX", LIS_rc%gridDesc(n,8))
         call PFIO_bundle%fmd_rst(n)%add_attribute("DY", LIS_rc%gridDesc(n,9))
      elseif (LIS_rc%lis_map_proj.eq."polar") then ! polar stereographic
         call PFIO_bundle%fmd_rst(n)%add_attribute("MAP_PROJECTION", "POLAR STEREOGRAPHIC")
         call PFIO_bundle%fmd_rst(n)%add_attribute("SOUTH_WEST_CORNER_LAT", SOUTH_WEST_CORNER_LAT)
         call PFIO_bundle%fmd_rst(n)%add_attribute("SOUTH_WEST_CORNER_LON", SOUTH_WEST_CORNER_LON)
         call PFIO_bundle%fmd_rst(n)%add_attribute("TRUELAT1", LIS_rc%gridDesc(n,10))
         call PFIO_bundle%fmd_rst(n)%add_attribute("ORIENT", LIS_rc%gridDesc(n,7))
         call PFIO_bundle%fmd_rst(n)%add_attribute("STANDARD_LON", LIS_rc%gridDesc(n,11))
         call PFIO_bundle%fmd_rst(n)%add_attribute("DX", LIS_rc%gridDesc(n,8))
         call PFIO_bundle%fmd_rst(n)%add_attribute("DY", LIS_rc%gridDesc(n,9))
      endif       

      ! Create the history collection
      PFIO_bundle%rst_id(n) = o_Clients%add_hist_collection(PFIO_bundle%fmd_rst(n))

   end subroutine PFIO_restart_metadata
!EOC
!------------------------------------------------------------------------------
   subroutine PFIO_restart_header(n, standard_name, long_name, units, vlevels, &
         valid_min, valid_max, var_flag)

! !USES: 

! !ARGUMENTS:     
      integer                    :: n
      character(len=*)           :: standard_name
      character(len=*)           :: long_name
      character(len=*)           :: units
      integer                    :: vlevels
      real                       :: valid_min
      real                       :: valid_max
      character(len=*), optional :: var_flag

      integer           :: deflate_level
      character(len=50) :: dim_names
      Type(Variable)    :: v
      character(len=50) :: var_flag_tmp


      deflate_level = LIS_rc%nc_deflate_lvl
      var_flag_tmp = ""
      if (present(var_flag)) var_flag_tmp = TRIM(var_flag)

      if (vlevels.gt.1) then
         dim_names = "ntiles"
         if (present(var_flag)) dim_names = 'ntiles'//','//TRIM(var_flag_tmp)
      else
         dim_names = "ntiles"
      endif

      v = Variable(type=PFIO_REAL32, dimensions=TRIM(dim_names), deflation=deflate_level)
      call v%add_attribute("units", trim(units))
      call v%add_attribute("standard_name", trim(standard_name))
      call v%add_attribute("long_name", trim(long_name))
      call v%add_attribute("scale_factor", 1.0)
      call v%add_attribute("add_offset", 0.0)
      call v%add_attribute("missing_value", LIS_rc%udef)
      call v%add_attribute("_FillValue", LIS_rc%udef)
      call v%add_attribute("vmin", valid_min)
      call v%add_attribute("vmax", valid_max)
      call PFIO_bundle%fmd_rst(n)%add_variable(trim(standard_name), v)

   end subroutine PFIO_restart_header
   !------------------------------------------------------------------------------
   subroutine PFIO_restart_write_real(n, file_name, var_data, var_name, n_levs)

! !ARGUMENTS:
      integer,          intent(in) :: n, n_levs
      real,             intent(in) :: var_data(:,:)
      character(len=*), intent(in) :: var_name
      character(len=*), intent(in) :: file_name
!
! !LOCAL VARIABLES:
      real, allocatable :: var_tile(:)
      integer           :: r, c, gid, tid, stid
      integer           :: peIdx, count1, ntiles


      !         call map_1dtile_to_1darray(n, var_data, var2d(:,1:n_levs), i1, i2, n_levs)
      !         ref =  ArrayReference(var2d(:,1:nlev))
      !         call o_Clients%collective_stage_data(PFIO_bundle%hist_id(n,vcol_id), &
      !                        TRIM(file_name), TRIM(var_name), ref, &
      !                        start        = [ i1, 1 ], &
      !                        global_start = [  1, 1 ], &
      !                        global_count = [ LIS_rc%glbntiles_red(n), nlev ] )


   end subroutine PFIO_restart_write_real
#endif
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
! -----------------------------------------------------------------------
   subroutine increment_field_counter(non_model_fields)
      integer, optional :: non_model_fields

      integer       :: nmodel_status
      nmodel_status = 0
      if (present(non_model_fields)) then
         nmodel_status = non_model_fields
      endif

      SELECT CASE(LIS_rc%wopt)
      CASE("1d tilespace")
         idx_field2d = idx_field2d + 1
      CASE("2d gridspace")
         idx_field3d = idx_field3d + 1
      CASE("2d ensemble gridspace")
         if (nmodel_status > 0) then
            idx_field3d = idx_field3d + 1
         else
            idx_field4d = idx_field4d + 1
         endif
      END SELECT
   end subroutine increment_field_counter
   ! -----------------------------------------------------------------------
   ! -----------------------------------------------------------------------
   !  subroutine set_field_index(local_var, idx_field3d, idx_field4d)
   !      integer, intent(in) :: idx_field3d, idx_field4d
   !      type(local_field_var), intent(inout) :: local_var(:)
   !
   !      SELECT CASE(LIS_rc%wopt)
   !      CASE("1d tilespace")
   !          CONTINUE
   !      CASE("2d gridspace")
   !          local_var(idx_field3d)%idx_field3d = idx_field3d
   !      CASE("2d ensemble gridspace")
   !          local_var(idx_field3d)%idx_field3d = idx_field3d
   !          local_var(idx_field4d)%idx_field4d = idx_field4d
   !      END SELECT
   !  end subroutine set_field_index
   ! -----------------------------------------------------------------------
   ! -----------------------------------------------------------------------
   subroutine allocate_local_var(n, PFIOmodel_idx, &
         local_var2D, local_var3D, local_var4D, i1, i2, j1, j2)
      integer, intent(in) :: n, PFIOmodel_idx
      integer, intent(in) :: i1, i2, j1, j2
      type(local_2Dfield_var), intent(inout) :: local_var2D(:)
      type(local_3Dfield_var), intent(inout) :: local_var3D(:)
      type(local_4Dfield_var), intent(inout) :: local_var4D(:)
      integer :: i

      SELECT CASE(LIS_rc%wopt)
      CASE("1d tilespace")
         DO i= 1, total_num_2Dfields(PFIOmodel_idx)
            allocate(local_var2D(i)%var2d(i1:i2,MAX_NUM_VLEVELS))
            local_var2D(i)%var2d = pfio_missing_value
         ENDDO
      CASE("2d gridspace")
         DO i= 1, total_num_3Dfields(PFIOmodel_idx)
            allocate(local_var3D(i)%var3d(i1:i2,j1:j2,MAX_NUM_VLEVELS))
            local_var3D(i)%var3d = pfio_missing_value
         ENDDO
      CASE("2d ensemble gridspace")
         DO i= 1, total_num_3Dfields(PFIOmodel_idx)
            allocate(local_var3D(i)%var3d(i1:i2,j1:j2,MAX_NUM_VLEVELS))
            local_var3D(i)%var3d = pfio_missing_value
         ENDDO
         DO i= 1, total_num_4Dfields(PFIOmodel_idx)
            allocate(local_var4D(i)%var4d(i1:i2,j1:j2,LIS_rc%nensem(n),MAX_NUM_VLEVELS))
            local_var4D(i)%var4d = pfio_missing_value
         ENDDO
      END SELECT
   end subroutine allocate_local_var
   ! -----------------------------------------------------------------------
   ! -----------------------------------------------------------------------
   subroutine deallocate_local_var(PFIOmodel_idx, local_var2D, local_var3D, local_var4D)
      integer,                 intent(in)    :: PFIOmodel_idx
      type(local_2Dfield_var), intent(inout) :: local_var2D(:)
      type(local_3Dfield_var), intent(inout) :: local_var3D(:)
      type(local_4Dfield_var), intent(inout) :: local_var4D(:)
      integer :: i

      SELECT CASE(LIS_rc%wopt)
      CASE("1d tilespace")
         DO i= 1, total_num_2Dfields(PFIOmodel_idx)
            deallocate(local_var2D(i)%var2d)
         ENDDO
         !deallocate(local_var2D)
      CASE("2d gridspace")
         DO i= 1, total_num_3Dfields(PFIOmodel_idx)
            deallocate(local_var3D(i)%var3d)
         ENDDO
         !deallocate(local_var3D)
      CASE("2d ensemble gridspace")
         DO i= 1, total_num_3Dfields(PFIOmodel_idx)
            deallocate(local_var3D(i)%var3d)
         ENDDO
         !deallocate(local_var3D)
         DO i= 1, total_num_4Dfields(PFIOmodel_idx)
            deallocate(local_var4D(i)%var4d)
         ENDDO
         !deallocate(local_var4D)
      END SELECT

      !DEALLOCATE(local_var)
   end subroutine deallocate_local_var
   ! -----------------------------------------------------------------------
   ! map_1dtile_to_2darray
   !    Map 1D local tile to a 2D local array
   ! -----------------------------------------------------------------------
   subroutine map_1dtile_to_2darray(n, var_tile, var3d, i1, i2, j1, j2, nlev)
      integer, intent(in) :: n
      integer, intent(in) :: i1, i2, j1, j2, nlev
      real, intent(in) :: var_tile(:,:)
      real, intent(out) :: var3d(i1:i2,j1:j2, nlev)

      integer :: r, c, t, peIdx
      integer :: ntiles, count1, gid

      var3d = pfio_missing_value ! LIS_rc%udef

      IF (LIS_rc%ntiles(n) .GT. 0) THEN
         peIdx = LIS_localPet+1
         count1 = 1
         do r=LIS_nss_halo_ind(n,peIdx),LIS_nse_halo_ind(n,peIdx)
            do c=LIS_ews_halo_ind(n,peIdx),LIS_ewe_halo_ind(n,peIdx)
               gid = c+(r-1)*LIS_rc%gnc(n)
               ntiles = LIS_domain(n)%ntiles_pergrid(gid)
               if (ntiles .ne. 0) then
                  IF (r.ge.LIS_nss_ind(n,peIdx) .and. r.le.LIS_nse_ind(n,peIdx) .and.&
                     c.ge.LIS_ews_ind(n,peIdx) .and. c.le.LIS_ewe_ind(n,peIdx)) THEN !points not in halo
                     var3d(c,r,1:nlev) = var_tile(count1,1:nlev)
                  END IF
                  count1 = count1 + 1
               endif
            enddo
         enddo
      ENDIF

   end subroutine map_1dtile_to_2darray
   ! -----------------------------------------------------------------------
   ! map_1dtile_to_1darray
   !    Map 1D local tile to a 2D local array
   ! -----------------------------------------------------------------------
   subroutine map_1dtile_to_1darray(n, var_tile, var2d, i1, i2, nlev)
      integer, intent(in) :: n
      integer, intent(in) :: i1, i2, nlev
      real, intent(in) :: var_tile(:,:)
      real, intent(out) :: var2d(i1:i2, nlev)

      integer :: r, c, t, l
      integer :: ntiles, count1, gid, stid, tid

      var2d = pfio_missing_value ! LIS_rc%udef

      l = LIS_localPet+1
      count1 = 1
      do r=LIS_nss_halo_ind(n,l),LIS_nse_halo_ind(n,l)
         do c=LIS_ews_halo_ind(n,l),LIS_ewe_halo_ind(n,l)
            gid = c+(r-1)*LIS_rc%gnc(n)
            ntiles = LIS_domain(n)%ntiles_pergrid(gid)
            stid = LIS_domain(n)%str_tind(gid)
            IF (r.ge.LIS_nss_ind(n,l) .and. r.le.LIS_nse_ind(n,l) .and.&
               c.ge.LIS_ews_ind(n,l) .and. c.le.LIS_ewe_ind(n,l)) THEN !points not in halo
               do t=1,ntiles
                  tid = stid + t-1
                  var2d(tid,1:nlev) = var_tile(count1,1:nlev)
                  count1 = count1 + 1
               enddo
            ELSE
               count1 = count1 + ntiles
            END IF
         enddo
      enddo

   end subroutine map_1dtile_to_1darray
   ! -----------------------------------------------------------------------
   ! get_LIS_interior_tile
   !    Determine the interior tiles for each domain
   ! -----------------------------------------------------------------------
   subroutine get_LIS_interior_tile(n, i1, i2)
      integer, intent(in) :: n
      integer, intent(out) :: i1, i2

      integer :: r, c, l, i1_n, i2_n, ntiles, stid
      integer :: count1, gid, tid, mymin, mymax,t
      logical :: first_time

      i1 = LIS_toffsets(n, LIS_localPet) + 1
      i2 = i1 + LIS_rc%ntiles(n) - 1

      first_time = .TRUE.
      l = LIS_localPet+1
      count1 = 0
      do r=LIS_nss_halo_ind(n,l),LIS_nse_halo_ind(n,l)
         do c=LIS_ews_halo_ind(n,l),LIS_ewe_halo_ind(n,l)
            gid = c+(r-1)*LIS_rc%gnc(n)
            IF (first_time) THEN
               i1_n = LIS_domain(n)%str_tind(gid)
               first_time = .FALSE.
            ENDIF
            count1 = count1 + LIS_domain(n)%ntiles_pergrid(gid)
         enddo
      enddo

      mymax = -9999
      mymin = 100000000
      count1 = 1
      do r=LIS_nss_halo_ind(n,l),LIS_nse_halo_ind(n,l)
         do c=LIS_ews_halo_ind(n,l),LIS_ewe_halo_ind(n,l)
            gid = c+(r-1)*LIS_rc%gnc(n)
            ntiles = LIS_domain(n)%ntiles_pergrid(gid)
            stid = LIS_domain(n)%str_tind(gid)
            IF (r.ge.LIS_nss_ind(n,l) .and. r.le.LIS_nse_ind(n,l) .and.&
               c.ge.LIS_ews_ind(n,l) .and. c.le.LIS_ewe_ind(n,l)) THEN !points not in halo
               do t=1,ntiles
                  tid = stid + t-1
                  if (tid .GE. mymax) mymax = tid
                  if (tid .LE. mymin) mymin = tid
                  count1 = count1 + 1
               enddo
            ELSE
               count1 = count1 + ntiles
            END IF
         enddo
      enddo
      i1 = i1_n
      i2 = i1 + count1 - 1

   end subroutine get_LIS_interior_tile
   !---------------------------------------------------------------------------
   ! Get the interior grid
   ! -----------------------------------------------------------------------
   subroutine get_LIS_interior_grid(n, i1, i2, j1, j2)
      integer, intent(in)  :: n
      integer, intent(out) :: i1, i2, j1, j2
      integer :: peIdx

      peIdx = LIS_localPet+1
      i1 = LIS_ews_ind(n,peIdx) ! lower longitude
      i2 = LIS_ewe_ind(n,peIdx) ! upper longitude
      j1 = LIS_nss_ind(n,peIdx) ! lower latitude
      j2 = LIS_nse_ind(n,peIdx) ! upper latitude

   end subroutine get_LIS_interior_grid
   ! -----------------------------------------------------------------------
   ! Get the global dimensions
   ! -----------------------------------------------------------------------
   subroutine get_LIS_globaldim(n, global_dim)
      integer, intent(in)  :: n
      integer, intent(out) :: global_dim(2)

      global_dim(1) = LIS_rc%gnc(n) ! longitudes / columns
      global_dim(2) = LIS_rc%gnr(n) ! latitudes / rows

   end subroutine get_LIS_globaldim
   ! -----------------------------------------------------------------------

#endif
end module LIS_PFIO_historyMod
