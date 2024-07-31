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
! !MODULE: climatology_forcingMod
!  \label{climatology_forcingMod}
!
! !REVISION HISTORY:
!  26 Sep 2016: KR Arsenault; Initial Specification 
!
! !NOTES: The implementation has not been tested for sub-hourly 
!         forcing climatology
! !INTERFACE:
module climatology_forcingMod
!
! !USES:
#if (defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif

  use LIS_constantsMod, only : LIS_CONST_PATH_LEN

  implicit none

!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: init_climatology      ! Defines the native resolution of 
                                  ! the input data
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  public :: clim_struc
!EOP

  type climatology_type_dec
     
     integer        :: nc
     integer        :: nr
     integer        :: ntimes
     real           :: ts
     real*8         :: findtime1, metforc_time1
     real*8         :: findtime2, metforc_time2

     character(len=LIS_CONST_PATH_LEN) :: directory
     character(50)  :: proj_name
     integer        :: proj_index
 
     logical        :: reset_flag
     logical        :: zterp_flags !(:) Allocatable later?
     
     real, allocatable :: metdata1(:,:) 
     real, allocatable :: metdata2(:,:) 

  end type climatology_type_dec
  
  type(climatology_type_dec) :: clim_struc

contains
  
!BOP
!
! !ROUTINE: init_climatology
!  \label{init_climatology}
! !INTERFACE:
  subroutine init_climatology(findex)

! !USES: 
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
    use netcdf
#endif
    use LIS_coreMod,    only : LIS_rc, LIS_domain
    use LIS_timeMgrMod, only : LIS_date2time, LIS_update_timestep,&
                               LIS_seconds2time
    use LIS_logMod,     only : LIS_logunit, LIS_verify, LIS_endrun
    use climatology_VariablesMod, only : climatology_Variables_init
    use climatology_SpatialInterpMod
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
    real, allocatable :: lat(:,:)
    real, allocatable :: lon(:,:)

    integer        :: numAtts
    character(40)  :: varname
    character(100) :: zterpflag_string

    real*8   :: LIS_time
    real     :: gmt
    integer  :: doy

! ______________________________________________________

    ! Forecast mode -- NOT Available at this time for this forcing reader:
    if( LIS_rc%forecastMode.eq.1 ) then
       write(LIS_logunit,*) '[ERR] Currently the climatology forcing reader'
       write(LIS_logunit,*) '[ERR]  is not set up to run in forecast mode.'
       write(LIS_logunit,*) '[ERR]  May be added in future releases.'
       write(LIS_logunit,*) '[ERR]  LIS forecast run-time ending.'
       call LIS_endrun()
    endif

    ! Initialize Climatology flag and times:
    gridDesci = 0.0
    clim_struc%reset_flag  = .false.
    clim_struc%metforc_time1 = 0.
    clim_struc%metforc_time2 = 0.

    ! Read lis.config entries:
    call readcrd_climatology()

    write(LIS_logunit,*) "[INFO] -- Obtain Climatology Forcing Dataset Parameters -- "

    ! Using starting LIS date/time values, locate starting climatology file: 
    fullfilename = "none"
    LIS_time = 0.
    call LIS_date2time( LIS_time, doy, gmt, &          ! Output
             LIS_rc%syr,LIS_rc%smo,LIS_rc%sda,0,0,0)   ! Input
    call get_climatology_filename( doy, &
             clim_struc%directory, fullfilename  )
    inquire( file=trim(fullfilename), exist=file_exists)
    if( file_exists ) then
       write(LIS_logunit,*) "[INFO] Parameters from forcing climatology file: ", &
             trim(fullfilename)
    else
       write(LIS_logunit,*) "[WARN] Missing Climatology forcing file: ", &
             trim(fullfilename)
       write(LIS_logunit,*) "[WARN]  LIS endrun being called."
       call LIS_endrun()
    endif

#if (defined USE_NETCDF3 || defined USE_NETCDF4)

!-- Read in initial Netcdf file .... 
    ios = nf90_open(path=fullfilename,&
          mode=NF90_NOWRITE,ncId=nid)
    call LIS_verify(ios,'Error in nf90_open in init_climatology')

  ! Obtain the time step:
    clim_struc%ts = 3600.    ! Temporary climatology timestep

#if 0
! LIS - finest timestep for runtime:  LIS_rc%nts(n)
!  To add climatology forcing timestep information get ...
    clim_struc%ts = 0.
    ios = nf90_inq_varid(nid,'time',timeId)
    call LIS_verify(ios,'time field not found in the metforcing climatology file')
    call LIS_verify(nf90_get_att(nid, timeID,&
        "time_increment", timeInc),&
        'nf90_get_att for time_increment failed in init_climatology')
    read(timeInc,*) clim_struc%ts
    call LIS_update_timestep(LIS_rc, 1, clim_struc%ts)
    write(LIS_logunit,*) "[INFO] LDT-Climatology Forcing Timestep: ",clim_struc%ts
#endif

  ! Number of hours per julian day ("ntimes"):
    clim_struc%ntimes = 0
    ios = nf90_inq_dimid(nid,"ntimes",ntimesId)
    call LIS_verify(ios,'Error in nf90_inq_dimid(ntimes) in init_climatology')
    ios = nf90_inquire_dimension(nid, ntimesId, len=clim_struc%ntimes)
    call LIS_verify(ios,'Error in nf90_inquire_dimension(ntimes) in init_climatology')
    write(LIS_logunit,*) "[INFO] LDT-Climatology Forcing Number of hours per day: ", &
         clim_struc%ntimes
   
    clim_struc%ts = 86400/clim_struc%ntimes

  ! Number of columns and rows:
    clim_struc%nc = 0
    clim_struc%nr = 0
    ios = nf90_inq_dimid(nid,"east_west",ncId)
    call LIS_verify(ios,'Error in nf90_inq_dimid(east_west) in init_climatology')
    ios = nf90_inq_dimid(nid,"north_south",nrId)
    call LIS_verify(ios,'Error in nf90_inq_dimid(north_south) in init_climatology')

    ios = nf90_inquire_dimension(nid, ncId, len=clim_struc%nc)  
    call LIS_verify(ios,'Error in nf90_inquire_dimension(nc) in init_climatology')
    ios = nf90_inquire_dimension(nid, nrId, len=clim_struc%nr)  
    call LIS_verify(ios,'Error in nf90_inquire_dimension(nr) in init_climatology')
    write(LIS_logunit,*) "[INFO] LDT-Climatology Forcing Number of Cols, Rows: ", &
                         clim_struc%nc, clim_struc%nr

  ! Obtain the lat and lon fields:
    allocate(lat(clim_struc%nc, clim_struc%nr))
    allocate(lon(clim_struc%nc, clim_struc%nr))

    ios = nf90_inq_varid(nid,'lat',latId)
    call LIS_verify(ios,'lat field not found in the metforcing climatology file')
    ios = nf90_inq_varid(nid,'lon',lonId)
    call LIS_verify(ios,'lon field not found in the metforcing climatology file')

    ios = nf90_get_var(nid,latId,lat)
    call LIS_verify(ios,'Error in nf90_get_var for lat in init_climatology')
    ios = nf90_get_var(nid,lonId,lon)
    call LIS_verify(ios,'Error in nf90_get_var for lon in init_climatology')

 !- Map Projection:
    ios = nf90_get_att(nid, NF90_GLOBAL, &
         'MAP_PROJECTION', clim_struc%proj_name )
    call LIS_verify(ios,'Error in nf90_get_att for MAP_PROJECTION in init_climatology')

  ! Define projection index and input GridDesc array:
    select case ( clim_struc%proj_name )

     case( "EQUIDISTANT CYLINDRICAL" )    

       write(LIS_logunit,*) "[INFO] LDT-Climatology Forcing Projection:  Latlon "
       clim_struc%proj_index = 0
       LIS_rc%met_proj(findex) = "latlon"
     
     ! Lower-left (ll) lat and long point:
       ios = nf90_get_att(nid, NF90_GLOBAL, &
            'SOUTH_WEST_CORNER_LAT', lllat )
       call LIS_verify(ios,'Error in nf90_get_att:SOUTH_WEST_CORNER_LAT in init_climatology')
       ios = nf90_get_att(nid, NF90_GLOBAL, &
            'SOUTH_WEST_CORNER_LON', lllon )
       call LIS_verify(ios,'Error in nf90_get_att:SOUTH_WEST_CORNER_LON in init_climatology')

       urlat = lat(clim_struc%nc, clim_struc%nr)
       urlon = lon(clim_struc%nc, clim_struc%nr)

     ! Resolutions (xres) and (yres):
       ios = nf90_get_att(nid, NF90_GLOBAL, 'DX', xres )
       call LIS_verify(ios,'Error in nf90_get_att:DX in init_climatology')
       ios = nf90_get_att(nid, NF90_GLOBAL, 'DY', yres )
       call LIS_verify(ios,'Error in nf90_get_att:DY in init_climatology')

       gridDesci = 0
       gridDesci(1)  = clim_struc%proj_index
       gridDesci(2)  = clim_struc%nc
       gridDesci(3)  = clim_struc%nr
       gridDesci(4)  = lllat
       gridDesci(5)  = lllon
       gridDesci(6)  = 128      
       gridDesci(7)  = urlat
       gridDesci(8)  = urlon
       gridDesci(9)  = xres
       gridDesci(10) = yres
       gridDesci(20) = 64     

     case default
        write(LIS_logunit,*) "[ERR] No other metforcing climatology projections "
        write(LIS_logunit,*) "[ERR]  are supported at this time ... only latlon."
        write(LIS_logunit,*) "[ERR] End run being called ..." 
        call LIS_endrun
    end select   

!- Extract zenith interp flags from global attribute:
!    ios = nf90_get_att(nid, NF90_GLOBAL, 'zenith_interp', zterpflag_string )
!    call LIS_verify(ios,'Error in nf90_get_att in init_climatology')
!!    zterp_flags = zterpflag_string ....  ! Still need to implement here ...
!    if( index(zterpflag_string,"false") > 0 ) then
!       clim_struc%zterp_flags = .false.
!    elseif( index(zterpflag_string,"true") > 0 ) then
!       clim_struc%zterp_flags = .true.
!    endif


!- Estimate number of forcing variables being readin from netcdf file:
    LIS_rc%met_nf(findex) = 0  

    call climatology_Variables_init( findex, nid, &
            clim_struc%nc, clim_struc%nr, &
            clim_struc%ntimes, LIS_rc%met_nf(findex)  )

    write(LIS_logunit,*) "[INFO] Number of metforcing climatology fields: ",&
          LIS_rc%met_nf(findex)

!- Close netCDF file.
   ios=nf90_close(nid)
   call LIS_verify(ios,'Error in nf90_close in init_climatology')
#endif


!- Allocate and calculate the interp points and weights to be
!   be used later when the metforcing data needs to be reprojected
!   onto the LIS run domain.

   write(LIS_logunit,*) "[INFO] -- NOTE :: "
   write(LIS_logunit,*) "[INFO]  Metforcing climatology reader at this time "
   write(LIS_logunit,*) "[INFO]  does not support time interpolation. To be "
   write(LIS_logunit,*) "[INFO]  implemented and released in the near future. "

   do n = 1, LIS_rc%nnest

       allocate(clim_struc%metdata1(LIS_rc%met_nf(findex),&
            LIS_rc%ngrid(n)))
       allocate(clim_struc%metdata2(LIS_rc%met_nf(findex),&
            LIS_rc%ngrid(n)))

       clim_struc%metdata1 = 0
       clim_struc%metdata2 = 0

      call climatology_interp_init( n, findex, gridDesci )
   end do
!--

  end subroutine init_climatology

end module climatology_forcingMod

