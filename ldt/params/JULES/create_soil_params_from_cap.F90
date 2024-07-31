!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LDT_misc.h"
!BOP
!
! !ROUTINE: create_soil_params_from_cap
! \label{create_soil_params_from_cap}
!
! !REVISION HISTORY:
!  03 Sept 2004: Sujay Kumar; Initial Specification
!  03  Oct 2013: KR Arsenault;  Expanded to read in native STATSGO v1 files
!  10  Apr 2017: Shugong Wang; Read CAP soil texture fraction data and calculate
!                soil hyraulic parameters for JULES model 
! !INTERFACE:
subroutine create_soil_params_from_cap( n, nc_cap_param, soil_type) 

! !USES:
  use LDT_coreMod, only : LDT_rc
  use LDT_logMod,  only : LDT_logunit,  LDT_endrun
  use LDT_gridmappingMod    
  use LDT_paramDataMod
  use LDT_fileIOMod

  implicit none
! !ARGUMENTS: 
  integer,intent(in)    :: n
  character(len=128)    :: nc_cap_param    ! NetCDF file of CAP soil texture fraction file 
  real, intent(inout)   :: soil_type(LDT_rc%lnc(n),LDT_rc%lnr(n))
!
! !DESCRIPTION:
!  This subroutine read the soil texture fraction data and derive soil hydraulic parameters.
!  The fraction data is in lat-lon coordinate and at 1/24 degree resolution. 
!
!  The arguments are:
!  \begin{description}
!  \item[n]
!   index of the nest
!  \item[nc_cap_param]
!   file name of the NetCDF CAP soil texture fraction data 
!  \end{description}
!EOP
   
   integer, parameter :: input_cols = 8640
   integer, parameter :: input_rows = 4320
   real,    parameter :: input_xres = 1.0/24.0
   real,    parameter :: input_yres = 1.0/24.0

   integer*4, parameter :: soil_nlayers = 1
 
   integer :: c, r, i, t 
   integer :: ftn
   integer :: length, nrec
   logical :: file_exists
   integer :: mi                       ! Total number of input param grid fgrd points
   integer :: mo                       ! Total number of output LIS grid fgrd points
   real    :: param_gridDesc(20)       ! Input parameter grid desc fgrd
   real,    allocatable  :: gi(:)      ! Input parameter 1d grid
   logical*1,allocatable :: li(:)      ! Input logical mask (to match gi)
   real      :: go1(LDT_rc%lnc(n)*LDT_rc%lnr(n))            ! Output lis 1d grid
   logical*1 :: lo1(LDT_rc%lnc(n)*LDT_rc%lnr(n))            ! Output logical mask (to match go)
   integer :: glpnc, glpnr               ! Parameter (global) total columns and rows
   integer :: subpnc, subpnr             ! Parameter subsetted columns and rows
   real    :: subparam_gridDesc(20)      ! Subsetted parameter grid desc array
   real, allocatable :: subset_usda(:,:)  ! Read input parameter
   integer, allocatable  :: lat_line(:,:), lon_line(:,:)


   real, allocatable:: silt_map(:,:), sand_map(:,:), clay_map(:,:)
   real, allocatable:: usda_type(:,:) 
! ___________________________________________________________________
   LDT_rc%waterclass = 0. 
   soil_type = LDT_rc%waterclass

!- Set parameter grid fgrd inputs:
!  Below grid based on info from:
!   http://dbwww.essc.psu.edu/dbtop/.link/1995-0795/proj_geo.html
   param_gridDesc(1)  = 0.                ! Latlon
   param_gridDesc(2)  = input_cols
   param_gridDesc(3)  = input_rows
   param_gridDesc(4)  =  -90.0+(1/48)   ! LL lat
   param_gridDesc(5)  = -180.0+(1/48)   ! LL lon *** domain starts from east 0 degree
   param_gridDesc(6)  = 128             ! ask Jim
   ! I believe that you should specify griddesc(4) = 0+(1/48) and griddesc(8) = 0 - (1/48).
   param_gridDesc(7)  =   90.0-(1/48)   ! UR lat
   param_gridDesc(8)  =  180.0-(1/48)   ! UR lon *** domain ends at west 0 degree 
   param_gridDesc(9)  = input_yres      ! dy: 0.0083333
   param_gridDesc(10) = input_xres      ! dx: 0.0083333
   param_gridDesc(20) = 64              ! indicating the starting lon is -180, double check 
! James V. Geiger
! griddesc(20) is used to specify reading columnwise first or rowwise first.  64 means columnwise first.  This is how most of our
! data are read in.
   LDT_rc%soiltext_gridDesc(n,:) = param_gridDesc(:)

   inquire(file=trim(nc_cap_param), exist=file_exists)
   if(.not.file_exists) then 
      write(LDT_logunit,*) "JULES soil fraction data ",trim(LDT_rc%txtfile(n))," not found."
      write(LDT_logunit,*) "Program stopping ..."
      call LDT_endrun
   endif

   write(unit=LDT_logunit,fmt=*) "[INFO] Reading JULES soil fraction file: ",&
        trim(nc_cap_param)

   call LDT_checkDomainExtents(n, param_gridDesc(:))

   allocate(silt_map(input_cols, input_rows)) 
   allocate(sand_map(input_cols, input_rows)) 
   allocate(clay_map(input_cols, input_rows)) 
   allocate(usda_type(input_cols, input_rows)) 

!- Read-in the JULSE CAP sand, silt and clay fraction data 
   call read_cap_soil_data(nc_cap_param, silt_map, sand_map, clay_map)

!  Create USDA soil type based on CAP sand, silt, and clay fraction data 
   call map_usda_texture_type(input_cols, input_rows, sand_map, silt_map, clay_map, usda_type)
   !open(unit=1001, file='usda_fine.txt', form='formatted', status='new', action='write' )
   !do c=1, input_cols
   ! do r=1, input_rows
   !   write(unit=1001, fmt='(F15.2)', advance='no') usda_type(c, r)
   ! end do
   ! write(unit=1001, fmt='(A1)') ''
   !end do
   !close(unit=1001)
! -------------------------------------------------------------------
!    PREPARE SUBSETTED PARAMETER GRID FOR READING IN NEEDED DATA
! -------------------------------------------------------------------

!- Map Parameter Grid Info to LIS Target Grid/Projection Info -- 
   subparam_gridDesc = 0.
   call LDT_RunDomainPts( n, LDT_rc%soiltext_proj, param_gridDesc(:), &
            glpnc, glpnr, subpnc, subpnr, subparam_gridDesc, lat_line, lon_line )
   
   allocate( subset_usda(subpnc, subpnr) )
   subset_usda = LDT_rc%waterclass
!- Subset parameter read-in array:
   do r = 1, subpnr
      do c = 1, subpnc
         ! reverse y-axis 
         subset_usda(c,r) = usda_type(lon_line(c,r), glpnr-lat_line(c,r)+1)
      enddo
   enddo

! -------------------------------------------------------------------
!     AGGREGATING FINE-SCALE GRIDS TO COARSER LIS OUTPUT GRID
! -------------------------------------------------------------------
   mi = input_cols*input_rows
   mo = LDT_rc%lnc(n)*LDT_rc%lnr(n)
   
   allocate( li(mi), gi(mi) )
   ! need ot couble check the initialization, may not matter for this case 
   gi  = float(LDT_rc%waterclass)
   li  = .false.
   lo1 = .false.

!- Assign 2-D input soil text to 1-D for aggregation routines:
   ! 12-category usda types, no water type 
   ! there are 0 values fractions, set soil type -9999
   i = 0
   do r = 1, subpnr
      do c = 1, subpnc
         i = i + 1
         gi(i) = subset_usda(c,r)
         if( gi(i) .ne. LDT_rc%udef ) li(i) = .true.
      enddo
   enddo
   ! need to do 


!- Transform parameter from original grid to LIS output grid:
   LDT_rc%soiltext_gridtransform(n) = "mode"
   call LDT_transform_paramgrid(n, LDT_rc%soiltext_gridtransform(n), &
                      !param_gridDesc, mi, 1, gi, li, mo, go1, lo1 )
                      subparam_gridDesc, mi, 1, gi, li, mo, go1, lo1 )

!- Convert 1D count to 2D grid fgrds:
   i = 0
   do r = 1, LDT_rc%lnr(n)
      do c = 1, LDT_rc%lnc(n)
         i = i + 1
         soil_type(c,r) = go1(i)
      enddo
   enddo
   
!   open(unit=1001, file='usda_nldas.txt', form='formatted', status='new', action='write' )
!   do c=1, LDT_rc%lnr(n)
!    do r=1, LDT_rc%lnc(n)
!      write(unit=1001, fmt='(F15.2)', advance='no') soil_type(c, r)
!    end do
!    write(unit=1001, fmt='(A1)') ''
!   end do
!   close(unit=1001)

   deallocate( gi, li )


   deallocate(silt_map) 
   deallocate(sand_map) 
   deallocate(clay_map) 

   write(unit=LDT_logunit,fmt=*) "[INFO] Done creating JULES soil hydraulic parameters based on sand, silt and clay fraction data"
 end subroutine create_soil_params_from_cap

 subroutine read_cap_soil_data(nc_cap_param, silt_map, sand_map, clay_map)
   use netcdf
   implicit none 
   character(len=*) :: nc_cap_param
   real, intent(inout) :: silt_map(8640,4320), clay_map(8640,4320), sand_map(8640,4320)
   real, allocatable :: data_tmp(:,:) 
   integer :: ncid, varid
   allocate(data_tmp(8640, 4320)) 
   ! open the file, NF90_NOWRITE tells netCDF we want read-only access to the file
   call jules_check_nc(nf90_open( trim(nc_cap_param), NF90_NOWRITE, ncid))
 
   ! get the varid of the data variable for silt
   call jules_check_nc(nf90_inq_varid(ncid, "volume_fraction_of_silt_in_soil", varid))
   ! read the volumetric fraction of silt
   !call jules_check_nc(nf90_get_var(ncid, varid, silt_map))
   call jules_check_nc(nf90_get_var(ncid, varid, data_tmp))
   ! shift data to -180 longitude 
   silt_map(4321:8640,:)  = data_tmp(1:4320,:)
   silt_map(1:4320,:)     = data_tmp(4321:8640,:)
   
   ! get the varid of the data variable for sand
   call jules_check_nc(nf90_inq_varid(ncid, "volume_fraction_of_sand_in_soil", varid))
   ! read the volumetric fraction of sand
   call jules_check_nc(nf90_get_var(ncid, varid, data_tmp))
   sand_map(4321:8640,:)  = data_tmp(1:4320,:)
   sand_map(1:4320,:)     = data_tmp(4321:8640,:)
 
   ! get the varid of the data variable for clay
   call jules_check_nc(nf90_inq_varid(ncid, "volume_fraction_of_clay_in_soil", varid))
   ! read the volumetric fraction of clay
   call jules_check_nc(nf90_get_var(ncid, varid, data_tmp))
   clay_map(4321:8640,:)  = data_tmp(1:4320,:)
   clay_map(1:4320,:)     = data_tmp(4321:8640,:)
   
   call jules_check_nc(nf90_close(ncid))
   deallocate(data_tmp) 
 end subroutine read_cap_soil_data
 
 subroutine jules_check_nc(status)
   use netcdf 
   implicit none
   integer, intent (in) :: status

   if(status /= nf90_noerr) then
     print *, trim(nf90_strerror(status))
     stop "Stopped"
   end if
 end subroutine jules_check_nc
  
! USDA 12-category soil types: https://prod.nrcs.usda.gov/Internet/FSE_DOCUMENTS/nrcs142p2_053196.xls  
! 1  SAND            : ((silt + 1.5*clay) < 15)
! 2  LOAMY SAND      : ((silt + 1.5*clay >= 15) && (silt + 2*clay < 30))
! 3  SANDY LOAM      : ((clay >= 7 && clay < 20) && (sand > 52) && ((silt + 2*clay) >= 30) || (clay < 7 && silt < 50 && (silt+2*clay)>=30))
! 4  LOAM            : ((clay >= 7 && clay < 27) && (silt >= 28 && silt < 50) && (sand <= 52))
! 5  SILT LOAM       : ((silt >= 50 && (clay >= 12 && clay < 27)) || ((silt >= 50 && silt < 80) && clay < 12))
! 6  SANDY CLAY LOAM : ((clay >= 20 && clay < 35) && (silt < 28) && (sand > 45)) 
! 7  CLAY LOAM       : ((clay >= 27 && clay < 40) && (sand > 20 && sand <= 45))
! 8  SILTY CLAY LOAM : ((clay >= 27 && clay < 40) && (sand  <= 20))
! 9  SANDY CLAY      : (clay >= 35 && sand > 45)
! 10 SILTY CLAY      : (clay >= 40 && silt >= 40)
! 11 CLAY            : (clay >= 40 && sand <= 45 && silt < 40)
! 12 SILT            : (silt >= 80 && clay < 12)

subroutine map_usda_texture_type(nlon, nlat, sand_map, silt_map, clay_map, usda_type)
  implicit none 
  integer, intent(in) :: nlon, nlat
  real, intent(in) :: silt_map(8640,4320), clay_map(8640,4320), sand_map(8640,4320)
  real, intent(inout) :: usda_type(8640, 4320); 
  real :: silt, clay, sand, type_id
  integer :: r, c

  do c=1, nlon
    do r=1, nlat
      silt = silt_map(c,r)*100
      sand = sand_map(c,r)*100
      clay = clay_map(c,r)*100
      if((silt <0) .or. (sand<0) .or. (clay <0)) then 
         type_id = -9999
      ! 1  SAND            : ((silt + 1.5*clay) < 15)
      else if((silt + 1.5*clay) < 15) then
        type_id = 1
      ! 2  LOAMY SAND      : ((silt + 1.5*clay >= 15) && (silt + 2*clay < 30))
      else if((silt + 1.5*clay >= 15) .and. (silt + 2*clay < 30)) then
        type_id = 2
      ! 3  SANDY LOAM      : ((clay >= 7 && clay < 20) && (sand > 52) && ((silt + 2*clay) >= 30) || (clay < 7 && silt < 50 && (silt+2*clay)>=30))
      else if((clay >= 7 .and. clay < 20) .and. (sand > 52) .and. ((silt + 2*clay) >= 30) .or. (clay < 7 .and. silt < 50 .and. (silt+2*clay)>=30)) then
        type_id = 3
      ! 4  LOAM            : ((clay >= 7 && clay < 27) && (silt >= 28 && silt < 50) && (sand <= 52))
      else if((clay >= 7 .and. clay < 27) .and. (silt >= 28 .and. silt < 50) .and. (sand <= 52)) then
        type_id = 4
      ! 5  SILT LOAM       : ((silt >= 50 && (clay >= 12 && clay < 27)) || ((silt >= 50 && silt < 80) && clay < 12))
      else if((silt >= 50 .and. (clay >= 12 .and. clay < 27)) .or. ((silt >= 50 .and. silt < 80) .and. clay < 12)) then
        type_id = 5
      ! 6  SANDY CLAY LOAM : ((clay >= 20 && clay < 35) && (silt < 28) && (sand > 45)) 
      else if((clay >= 20 .and. clay < 35) .and. (silt < 28) .and. (sand > 45)) then
        type_id = 6
      ! 7  CLAY LOAM       : ((clay >= 27 && clay < 40) && (sand > 20 && sand <= 45))
      else if((clay >= 27 .and. clay < 40) .and. (sand > 20 .and. sand <= 45)) then
        type_id = 7
      ! 8  SILTY CLAY LOAM : ((clay >= 27 && clay < 40) && (sand  <= 20))
      else if((clay >= 27 .and. clay < 40) .and. (sand  <= 20)) then
        type_id = 8
      ! 9  SANDY CLAY      : (clay >= 35 && sand > 45)
      else if(clay >= 35 .and. sand > 45) then
        type_id = 9
      ! 10 SILTY CLAY      : (clay >= 40 && silt >= 40)
      else if(clay >= 40 .and. silt >= 40) then
        type_id = 10
      ! 11 CLAY            : (clay >= 40 && sand <= 45 && silt < 40)
      else if(clay >= 40 .and. sand <= 45 .and. silt < 40) then
        type_id = 11
      ! 12 SILT            : (silt >= 80 && clay < 12)
      else if (silt >= 80 .and. clay < 12) then
        type_id = 12
      else 
        write(*,*) "Error: uncategorized soil types: ", sand, silt, clay
      endif
      usda_type(c,r) = type_id
    enddo
  enddo
end subroutine map_usda_texture_type
