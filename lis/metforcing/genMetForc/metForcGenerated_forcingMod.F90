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
! !MODULE: metForcGenerated_forcingMod
!  \label{metForcGenerated_forcingMod}
!
! !REVISION HISTORY:
!  06 Jan 2015: KR Arsenault; Initial Specification 
!
! !INTERFACE:
module metForcGenerated_forcingMod
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
  public :: init_metForcGenerated      ! Defines the native resolution of 
                                       ! the input data
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  public :: metForcGen_struc
!EOP

  type metForcGenerated_type_dec

     integer        :: nc
     integer        :: nr
     real           :: ts
     real*8         :: metforc_time1
     real*8         :: metforc_time2

     character(len=LIS_CONST_PATH_LEN) :: directory
     character(50)  :: proj_name
     integer        :: proj_index
 
     logical        :: zterp_flags !(:) Allocatable later?

     integer           :: nIter, st_iterid,en_iterid

     real, allocatable :: metdata1(:,:,:) 
     real, allocatable :: metdata2(:,:,:)      

  end type metForcGenerated_type_dec
  
  type(metForcGenerated_type_dec) :: metForcGen_struc

contains
  
!BOP
!
! !ROUTINE: init_metForcGenerated
!  \label{init_metForcGenerated}
! !INTERFACE:
  subroutine init_metForcGenerated(findex)

! !USES: 
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
    use netcdf
#endif
    use LIS_coreMod,    only : LIS_rc, LIS_domain
    use LIS_timeMgrMod, only : LIS_date2time, LIS_update_timestep,&
                               LIS_seconds2time
    use LIS_logMod,     only : LIS_logunit, LIS_verify, LIS_endrun
    use metForcGen_VariablesMod, only : metForcGen_Variables_init
    use metForcGen_SpatialInterpMod
    use LIS_forecastMod
!EOP

    implicit none
    integer, intent(in) :: findex

    integer  :: n, i
    integer  :: ios
    integer  :: nid, ncId, nrId
    integer  :: timeId, lonId, latId
    integer  :: varid
    real     :: gridDesci(50)
    logical  :: file_exists
    character(len=LIS_CONST_PATH_LEN) :: fullfilename

    integer  :: seconds, da, hr, mn, ss, icount
    character*50 :: timeInc

    real     :: lllat, lllon
    real     :: urlat, urlon
    real     :: xres, yres
    real, allocatable :: lat(:,:)
    real, allocatable :: lon(:,:)

    integer        :: numAtts
    character(40)  :: varname
    character(100) :: zterpflag_string

! ______________________________________________________

    gridDesci = 0.0

    call readcrd_metForcGenerated()

    write(LIS_logunit,*) "[INFO] -- Obtain Forcing Dataset Parameters -- "

! - Need to loop over starting files to narrow down the correct file ...
    icount  = 0
    seconds = 0
    fullfilename = ""
    do
       icount  = icount + 1
!       seconds = seconds + 600  ! Advance every 10 minutes until file found
       seconds = seconds + 300  ! Advance every 5 minutes until file found
       call LIS_seconds2time( seconds, da, hr, mn, ss )
       call get_metForcGen_filename( 1, 1, findex, LIS_rc%yr, LIS_rc%mo, &
!            da, hr, 0, metForcGen_struc%directory, fullfilename  )
            LIS_rc%da, hr, 0, metForcGen_struc%directory, fullfilename  )
       inquire( file=trim(fullfilename), exist=file_exists)
       if( file_exists ) then
          write(LIS_logunit,*) "[INFO] -- Parameters from forcing file: ", &
               trim(fullfilename)
          exit
       endif
       if( icount > 1000 ) then
          write(LIS_logunit,*)"[ERR] No LDT-generated forcing file, "
          write(LIS_logunit,*) trim(fullfilename)
          write(LIS_logunit,*)"[ERR] is missing and can't extract file grid and time info."
          write(LIS_logunit,*)"[ERR] Calling End Run ... "
          call LIS_endrun
       endif
    end do

#if (defined USE_NETCDF3 || defined USE_NETCDF4)

!-- Read in initial Netcdf file .... 
    ios = nf90_open(path=fullfilename,&
          mode=NF90_NOWRITE,ncId=nid)
    call LIS_verify(ios,'Error in nf90_open in init_metForcGenerated')

  ! Obtain the time step:
    metForcGen_struc%ts = 0.
    ios = nf90_inq_varid(nid,'time',timeId)
    call LIS_verify(ios,'time field not found in the LDT-generated forcing file')

    call LIS_verify(nf90_get_att(nid, timeID,&
        "time_increment", timeInc),&
        'nf90_get_att for time_increment failed in init_metForcGenerated')
    read(timeInc,*) metForcGen_struc%ts
    call LIS_update_timestep(LIS_rc, 1, metForcGen_struc%ts)
    write(LIS_logunit,*) "[INFO] LDT-Generated Forcing Timestep: ",metForcGen_struc%ts

  ! Number of columns and rows:
    metForcGen_struc%nc = 0
    metForcGen_struc%nr = 0

    ios = nf90_inq_dimid(nid,"east_west",ncId)
    call LIS_verify(ios,'Error in nf90_inq_dimid in init_metForcGenerated')
    ios = nf90_inq_dimid(nid,"north_south",nrId)
    call LIS_verify(ios,'Error in nf90_inq_dimid in init_metForcGenerated')

    ios = nf90_inquire_dimension(nid, ncId, len=metForcGen_struc%nc)  
    call LIS_verify(ios,'Error in nf90_inquire_dimension in init_metForcGenerated')
    ios = nf90_inquire_dimension(nid, nrId, len=metForcGen_struc%nr)  
    call LIS_verify(ios,'Error in nf90_inquire_dimension in init_metForcGenerated')

    write(LIS_logunit,*) "[INFO] LDT-Generated Forcing Number of Cols, Rows: ", &
                         metForcGen_struc%nc, metForcGen_struc%nr

  ! Obtain the lat and lon fields:
    allocate(lat(metForcGen_struc%nc, metForcGen_struc%nr))
    allocate(lon(metForcGen_struc%nc, metForcGen_struc%nr))

    ios = nf90_inq_varid(nid,'lat',latId)
    call LIS_verify(ios,'lat field not found in the LDT-generated forcing file')
    ios = nf90_inq_varid(nid,'lon',lonId)
    call LIS_verify(ios,'lon field not found in the LDT-generated forcing file')

    ios = nf90_get_var(nid,latId,lat)
    call LIS_verify(ios,'Error in nf90_get_var for lat in init_metForcGenerated')
    ios = nf90_get_var(nid,lonId,lon)
    call LIS_verify(ios,'Error in nf90_get_var for lon in init_metForcGenerated')

 !- Map Projection:
    ios = nf90_get_att(nid, NF90_GLOBAL, &
         'MAP_PROJECTION', metForcGen_struc%proj_name )
    call LIS_verify(ios,'Error in nf90_get_att for MAP_PROJECTION in init_metForcGenerated')

  ! Define projection index and input GridDesc array:
    select case ( metForcGen_struc%proj_name )

     case( "EQUIDISTANT CYLINDRICAL" )    

       write(LIS_logunit,*) "[INFO] LDT-Generated Forcing Projection:  Latlon "
       metForcGen_struc%proj_index = 0
       LIS_rc%met_proj(findex) = "latlon"
     
     ! Lower-left (ll) lat and long point:
       ios = nf90_get_att(nid, NF90_GLOBAL, &
            'SOUTH_WEST_CORNER_LAT', lllat )
       call LIS_verify(ios,'Error in nf90_get_att:SOUTH_WEST_CORNER_LAT in init_metForcGenerated')
       ios = nf90_get_att(nid, NF90_GLOBAL, &
            'SOUTH_WEST_CORNER_LON', lllon )
       call LIS_verify(ios,'Error in nf90_get_att:SOUTH_WEST_CORNER_LON in init_metForcGenerated')

       urlat = lat(metForcGen_struc%nc, metForcGen_struc%nr)
       urlon = lon(metForcGen_struc%nc, metForcGen_struc%nr)

     ! Resolutions (xres) and (yres):
       ios = nf90_get_att(nid, NF90_GLOBAL, 'DX', xres )
       call LIS_verify(ios,'Error in nf90_get_att:DX in init_metForcGenerated')
       ios = nf90_get_att(nid, NF90_GLOBAL, 'DY', yres )
       call LIS_verify(ios,'Error in nf90_get_att:DY in init_metForcGenerated')

       gridDesci = 0
       gridDesci(1)  = metForcGen_struc%proj_index
       gridDesci(2)  = metForcGen_struc%nc
       gridDesci(3)  = metForcGen_struc%nr
       gridDesci(4)  = lllat
       gridDesci(5)  = lllon
       gridDesci(6)  = 128      
       gridDesci(7)  = urlat
       gridDesci(8)  = urlon
       gridDesci(9)  = xres
       gridDesci(10) = yres
       gridDesci(20) = 64     

     case default
        write(LIS_logunit,*) "[ERR] No other LDT-generated forcing projections "
        write(LIS_logunit,*) "[ERR]  are supported at this time ... only latlon."
        write(LIS_logunit,*) "[ERR] End run being called ..." 
        call LIS_endrun
    end select   

!- Extract zenith interp flags from global attribute:
    ios = nf90_get_att(nid, NF90_GLOBAL, 'zenith_interp', zterpflag_string )
    call LIS_verify(ios,'Error in nf90_get_att in init_metForcGenerated')
!    zterp_flags = zterpflag_string ....  ! Still need to implement here ...
    if( index(zterpflag_string,"false") > 0 ) then
       metForcGen_struc%zterp_flags = .false.
    elseif( index(zterpflag_string,"true") > 0 ) then
       metForcGen_struc%zterp_flags = .true.
    endif

!- Estimate number of forcing variables being readin from netcdf file:
    LIS_rc%met_nf(findex) = 0  

    call metForcGen_Variables_init( findex, nid, &
            metForcGen_struc%nc, metForcGen_struc%nr, &
            LIS_rc%met_nf(findex)  )

    write(LIS_logunit,*) "[INFO] Number of LDT-generated forcing fields: ",&
          LIS_rc%met_nf(findex)

!- Close netCDF file.
   ios=nf90_close(nid)
   call LIS_verify(ios,'Error in nf90_close in init_metForcGenerated')

#endif

!- Allocate and calculate the interp points and weights to be
!   be used later when the metforcing data needs to be reprojected
!   onto the LIS run domain.
   do n = 1, LIS_rc%nnest
      
      if(LIS_rc%forecastMode.eq.1) then

        if(mod(LIS_rc%nensem(n),&
             LIS_forecast_struc(1)%niterations).ne.0) then
           write(LIS_logunit,*) '[ERR] The number of ensembles must be a multiple'
           write(LIS_logunit,*) '[ERR] of the number of iterations '
           write(LIS_logunit,*) '[ERR] nensem = ',LIS_rc%nensem(n)
           write(LIS_logunit,*) '[ERR] niter = ',LIS_forecast_struc(1)%niterations
           call LIS_endrun()
        endif

        metForcGen_struc%st_iterid = LIS_forecast_struc(1)%st_iterId
        metForcGen_struc%en_iterId = LIS_forecast_struc(1)%niterations
        metForcGen_struc%nIter = LIS_forecast_struc(1)%niterations

        allocate(metForcGen_struc%metdata1(LIS_forecast_struc(1)%niterations,&
             LIS_rc%met_nf(findex),&
             LIS_rc%ngrid(n)))
        allocate(metForcGen_struc%metdata2(LIS_forecast_struc(1)%niterations,&
             LIS_rc%met_nf(findex),&
             LIS_rc%ngrid(n)))

      else   ! Non-forecast mode

        metForcGen_struc%st_iterid = 1
        metForcGen_struc%en_iterId = 1
        metForcGen_struc%nIter = 1
        
        allocate(metForcGen_struc%metdata1(1,&
             LIS_rc%met_nf(findex),&
             LIS_rc%ngrid(n)))
        allocate(metForcGen_struc%metdata2(1,&
             LIS_rc%met_nf(findex),&
             LIS_rc%ngrid(n)))
      endif

      metForcGen_struc%metdata1 = 0
      metForcGen_struc%metdata2 = 0

      call metForcGen_interp_init( n, findex, gridDesci )
   end do
!--

  end subroutine init_metForcGenerated

end module metForcGenerated_forcingMod

