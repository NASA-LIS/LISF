!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LIS_misc.h"
!BOP
!
! !MODULE: genEnsFcst_forcingMod
!  \label{genEnsFcst_forcingMod}
!
! !REVISION HISTORY:
!  16 Dec 2016: KR Arsenault; Initial Specification 
!
! !NOTES: The implementation has not been tested for sub-hourly 
!         forcing genEnsFcst
! !INTERFACE:
module genEnsFcst_forcingMod
!
! !USES:
#if (defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif
  
  use LIS_constantsMod, only : LIS_CONST_PATH_LEN
!
  implicit none

!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: init_genEnsFcst      ! Defines the native resolution of 
                                  ! the input data
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  public :: genensfcst_struc
!EOP

  type genEnsFcst_type_dec
     
     integer        :: nc               ! Number of columns (long points)
     integer        :: nr               ! Number of rows (lat points)
     integer        :: start_nc, start_nr  ! Starting nc / nr position for subset region
     integer        :: ntimes           ! Number of time points in a file
     integer        :: time_intv        ! Daily time interval, number of times in day
     real           :: ts               ! Forcing timestep
     real*8         :: findtime1, metforc_time1   ! File 1 flag and time (LIS)
     real*8         :: findtime2, metforc_time2   ! File 2 flag and time (LIS)

     character(len=LIS_CONST_PATH_LEN) :: directory        ! Directory path of where files reside
     character(50)  :: proj_name        ! Projection name
     integer        :: proj_index       ! Projection type index
 
     logical        :: reset_flag       ! Reset parameters flag
     logical        :: zterp_flags      !(:) Allocatable later?

     character(20)  :: fcst_type        ! Forecast type name
     character(40)  :: user_spec        ! Forecast file, user-specified naming convention
     integer        :: max_ens_members  ! Max number of forecast ensemble members
   
     ! Specify forecast initialization date:
     integer        :: fcst_inityr      ! Forecast initial year
     integer        :: fcst_initmo      ! Forecast initial month
     integer        :: fcst_initda      ! Forecast initial day

     real, allocatable :: metdata1(:,:,:)  ! Metforcing data from file/bookend 1
     real, allocatable :: metdata2(:,:,:)  ! Metforcing data from file/bookend 2
     
  end type genEnsFcst_type_dec
  
  type(genEnsFcst_type_dec) :: genensfcst_struc

contains
  
!BOP
!
! !ROUTINE: init_genEnsFcst
!  \label{init_genEnsFcst}
! !INTERFACE:
  subroutine init_genEnsFcst(findex)

! !USES: 
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
    use netcdf
#endif
    use LIS_coreMod,    only : LIS_rc, LIS_domain
    use LIS_timeMgrMod, only : LIS_date2time, LIS_update_timestep,&
                               LIS_seconds2time
    use LIS_logMod,     only : LIS_logunit, LIS_verify, LIS_endrun
    use genEnsFcst_VariablesMod, only : genEnsFcst_Variables_init
    use genEnsFcst_SpatialInterpMod
    use LIS_gridmappingMod
!EOP

  implicit none
  integer, intent(in) :: findex

  integer  :: n, i
  integer  :: ios
  integer  :: nid, ncId, nrId, ntimesId
  integer  :: timeId, lonId, latId
  integer  :: varid
  real     :: gridDesci(50)
  logical  :: file_exists
  character(len=LIS_CONST_PATH_LEN) :: fullfilename

  integer  :: da, hr, mn, ss
  character*50 :: timeInc
  real     :: lllat, lllon
  real     :: urlat, urlon
  real     :: xres, yres
  real, allocatable :: lat(:)
  real, allocatable :: lon(:)

  integer        :: numAtts
  character(40)  :: varname
  character(100) :: zterpflag_string

  integer :: tdel, tintv_validhr, hindex, tindex
  integer, dimension(12) :: mon_numdays
  data mon_numdays /31,28,31,30,31,30,31,31,30,31,30,31/

  integer :: subpnc, subpnr, glpnc, glpnr
  real    :: sub_gridDesci(20)
  integer, allocatable  :: lat_line(:,:), lon_line(:,:)
! _________________________________________________________

   ! Forecast mode -- NOT Available at this time for this forcing reader:
   if( LIS_rc%forecastMode.eq.1 ) then
      write(LIS_logunit,*)'[ERR] Currently the generic ensemble forecast reader'
      write(LIS_logunit,*)'[ERR] is not set up to run in an ESP-type forecast run mode.'
      write(LIS_logunit,*)'[ERR] May be added in future releases.'
      write(LIS_logunit,*)'[ERR] LIS-ESP forecast run-time ending.'
      call LIS_endrun()
   endif

   ! Leap year
   if( (mod(LIS_rc%yr,4).eq.0.and.mod(LIS_rc%yr,100).ne.0) .or.&
       (mod(LIS_rc%yr,400).eq.0) ) then
      mon_numdays(2) = 29
   else
      mon_numdays(2) = 28
   endif

   ! Initialize ensemble forecast dataset flag and times:
   gridDesci = 0.0
   genensfcst_struc%reset_flag  = .false.
   genensfcst_struc%metforc_time1 = 0.
   genensfcst_struc%metforc_time2 = 0.

   ! Read lis.config entries:
   call readcrd_genEnsFcst()

   write(LIS_logunit,*) "[INFO] -- Obtain ensemble forecast dataset parameters -- "
   write(LIS_logunit,*) "[INFO] -- Base forecast dataset name: ", &
      genensfcst_struc%fcst_type

   ! Locate starting genEnsFcst file: 
   fullfilename = "none"
   call get_genEnsFcst_filename( genensfcst_struc%fcst_type, &
      genensfcst_struc%user_spec, &
!      LIS_rc%syr, LIS_rc%smo, &    ! (Original code prior to 09-02-2022)
      genensfcst_struc%fcst_inityr, genensfcst_struc%fcst_initmo, &  ! New config entry
      1, LIS_rc%syr, LIS_rc%smo, &
      genensfcst_struc%directory, fullfilename )

   inquire( file=trim(fullfilename), exist=file_exists)
   if( file_exists ) then
      write(LIS_logunit,*) "[INFO] Parameters from forcing genEnsFcst file: ", &
            trim(fullfilename)
   else
      write(LIS_logunit,*) "[ERR] Missing ensemble forecast file: ", &
            trim(fullfilename)
      write(LIS_logunit,*) "[ERR] LIS endrun being called."
      call LIS_endrun()
   endif

!-- Read in initial Netcdf file .... 
#if (defined USE_NETCDF3 || defined USE_NETCDF4)
   ios = nf90_open(path=fullfilename,&
         mode=NF90_NOWRITE,ncId=nid)
   call LIS_verify(ios,'Error in nf90_open in init_genEnsFcst')

   ! Initialize forcing time step:
   genensfcst_struc%ts = 3600.    ! Temporary genEnsFcst timestep

   ! Number of time points in monthly forecast file ("ntimes"):
   genensfcst_struc%ntimes = 0
   ios = nf90_inq_dimid(nid,"time",ntimesId)
   call LIS_verify(ios, &
      'Error in nf90_inq_dimid(time) in init_genEnsFcst')
   ios = nf90_inquire_dimension(nid, ntimesId, len=genensfcst_struc%ntimes)
   call LIS_verify(ios, &
      'Error in nf90_inquire_dimension(time) in init_genEnsFcst')
   write(LIS_logunit,*) "[INFO] Ensemble forecast month number of timesteps: ",&
        genensfcst_struc%ntimes

   ! Determine number of time intervals in a day, starting by month, and
   genensfcst_struc%time_intv = genensfcst_struc%ntimes/mon_numdays(LIS_rc%smo)
   tdel = (24/genensfcst_struc%time_intv)
   ! Final LIS-ensemble forcing timestep:
   genensfcst_struc%ts = tdel*3600

   ! Determine final time index (tindex) of the point in the file to be read:
   hindex = (LIS_rc%hr/tdel)+1
   tindex = genensfcst_struc%time_intv * (LIS_rc%sda-1) + hindex
!   print *, "Tintv,tdel,hindex,tindex: ",genensfcst_struc%time_intv, tdel, hindex, tindex
 ! ___

   ! Number of columns and rows:
   genensfcst_struc%nc = 0
   genensfcst_struc%nr = 0
   ios = nf90_inq_dimid(nid,"longitude",ncId)
!   ios = nf90_inq_dimid(nid,"lon",ncId)
   call LIS_verify(ios, &
      'Error in nf90_inq_dimid(longitude_nc) in init_genEnsFcst')
   ios = nf90_inq_dimid(nid,"latitude",nrId)
!   ios = nf90_inq_dimid(nid,"lat",nrId)
   call LIS_verify(ios, &
      'Error in nf90_inq_dimid(latitude_nr) in init_genEnsFcst')

   ios = nf90_inquire_dimension(nid, ncId, len=genensfcst_struc%nc)  
   call LIS_verify(ios,'Error in nf90_inquire_dimension(nc) in init_genEnsFcst')
   ios = nf90_inquire_dimension(nid, nrId, len=genensfcst_struc%nr)  
   call LIS_verify(ios,'Error in nf90_inquire_dimension(nr) in init_genEnsFcst')
   write(LIS_logunit,*) &
      "[INFO] Ensemble forecast dataset Forcing Number of Cols, Rows: ", &
      genensfcst_struc%nc, genensfcst_struc%nr

   ! Obtain the lat and lon fields:
   allocate(lat(genensfcst_struc%nr))
   allocate(lon(genensfcst_struc%nc))

   ios = nf90_inq_varid(nid,'latitude',latId)
!   ios = nf90_inq_varid(nid,'lat',latId)
   call LIS_verify(ios, &
      'latitude var field not found in the metforcing genEnsFcst file')
   ios = nf90_inq_varid(nid,'longitude',lonId)
!   ios = nf90_inq_varid(nid,'lon',lonId)
   call LIS_verify(ios, &
      'longitude var field not found in the metforcing genEnsFcst file')

   ios = nf90_get_var(nid,latId,lat)
   call LIS_verify(ios, &
      'Error in nf90_get_var for latitude var in init_genEnsFcst')
   ios = nf90_get_var(nid,lonId,lon)
   call LIS_verify(ios, &
      'Error in nf90_get_var for longitude var in init_genEnsFcst')
 
   ! Map Projection:
   ios = nf90_get_att(nid, NF90_GLOBAL, &
         'MAP_PROJECTION', genensfcst_struc%proj_name )
   call LIS_verify(ios, &
      'Error in nf90_get_att for MAP_PROJECTION in init_genEnsFcst')

   ! Define projection index and input GridDesc array:
   select case ( genensfcst_struc%proj_name )

     case( "EQUIDISTANT CYLINDRICAL" )    

       write(LIS_logunit,*) &
          "[INFO] Ensemble forecast dataset Forcing Projection:  Latlon "
       genensfcst_struc%proj_index = 0
       LIS_rc%met_proj(findex) = "latlon"
     
       ! Lower-left (ll) lat and long point:
       ios = nf90_get_att(nid, NF90_GLOBAL, &
            'SOUTH_WEST_CORNER_LAT', lllat )
       call LIS_verify(ios, &
          'Error in nf90_get_att:SOUTH_WEST_CORNER_LAT in init_genEnsFcst')
       ios = nf90_get_att(nid, NF90_GLOBAL, &
            'SOUTH_WEST_CORNER_LON', lllon )
       call LIS_verify(ios, &
          'Error in nf90_get_att:SOUTH_WEST_CORNER_LON in init_genEnsFcst')

       ! Upper-right (ur) lat / long point:
       urlat = lat(genensfcst_struc%nr)
       urlon = lon(genensfcst_struc%nc)

       ! Resolutions (xres) and (yres):
       ios = nf90_get_att(nid, NF90_GLOBAL, 'DX', xres )
       call LIS_verify(ios,'Error in nf90_get_att:DX in init_genEnsFcst')
       ios = nf90_get_att(nid, NF90_GLOBAL, 'DY', yres )
       call LIS_verify(ios,'Error in nf90_get_att:DY in init_genEnsFcst')

       gridDesci = 0
       gridDesci(1)  = genensfcst_struc%proj_index
       gridDesci(2)  = genensfcst_struc%nc
       gridDesci(3)  = genensfcst_struc%nr
       gridDesci(4)  = lllat
       gridDesci(5)  = lllon
       gridDesci(6)  = 128      
       gridDesci(7)  = urlat
       gridDesci(8)  = urlon
       gridDesci(9)  = xres
       gridDesci(10) = yres
       gridDesci(20) = 64     

! -- IF WANT TO SUBSET THE INPUT NETCDF FILES IN FUTURE ....
! -- NEED TO ADD ADDITIONAL CODE IN LIS_RunDomainPts
! --  TO SUPPORT FORCING FIELDS COARSER THAN RUN DOMAIN ...
#if 0
       ! Number of columns and rows:
       sub_gridDesci = 0.
       subpnc = 0; subpnr = 0
       genensfcst_struc%start_nc = 0
       genensfcst_struc%start_nr = 0
       call LIS_RunDomainPts( 1, LIS_rc%met_proj(findex), gridDesci(:), &
                  glpnc, glpnr, subpnc, subpnr, &
                  sub_gridDesci, lat_line, lon_line )

       genensfcst_struc%start_nc = lon_line(1,1)
       genensfcst_struc%start_nr = lat_line(1,1)
       ! Subsetted domains:
       genensfcst_struc%nc = subpnc
       genensfcst_struc%nr = subpnr
#endif 

     case default
        write(LIS_logunit,*) "[ERR] No other metforcing genEnsFcst projections"
        write(LIS_logunit,*) "[ERR] are supported at this time ... only latlon."
        write(LIS_logunit,*) "[ERR] End run being called ..." 
        call LIS_endrun
    end select   

!- Estimate number of forcing variables being read in from netcdf file:
    LIS_rc%met_nf(findex) = 0  

    call genEnsFcst_Variables_init( findex, nid, &
            genensfcst_struc%nc, genensfcst_struc%nr, &
            1, LIS_rc%met_nf(findex) )    ! 1 - starting time index

    write(LIS_logunit,*) "[INFO] Number of Ensemble forecast fields: ",&
          LIS_rc%met_nf(findex)

   ! Close netCDF file.
   ios=nf90_close(nid)
   call LIS_verify(ios,'Error in nf90_close in init_genEnsFcst')
#endif

!- Allocate and calculate the interp points and weights to be
!   be used later when the metforcing data needs to be reprojected
!   onto the LIS run domain.

   write(LIS_logunit,*) "[INFO] -- NOTE ::"
   write(LIS_logunit,*) "[INFO] Metforcing 'genEnsFcst' reader at this time"
   write(LIS_logunit,*) "[INFO] does not support time interpolation. To be"
   write(LIS_logunit,*) "[INFO] implemented and released in the near future."

   LIS_rc%met_nensem(findex) = genensfcst_struc%max_ens_members

   do n = 1, LIS_rc%nnest
      call genEnsFcst_interp_init( n, findex, gridDesci )

      allocate(genensfcst_struc%metdata1(LIS_rc%met_nf(findex),&
               genensfcst_struc%max_ens_members,LIS_rc%ngrid(n)))
      allocate(genensfcst_struc%metdata2(LIS_rc%met_nf(findex),&
               genensfcst_struc%max_ens_members,LIS_rc%ngrid(n)))
   end do

   genensfcst_struc%metdata1 = 0.
   genensfcst_struc%metdata2 = 0.


  end subroutine init_genEnsFcst

end module genEnsFcst_forcingMod

