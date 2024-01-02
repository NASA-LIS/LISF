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
! !ROUTINE: create_soil_albedo_from_whs
! \label{create_soil_albedo_from_whs}
!
! !REVISION HISTORY:
!  03 Sept 2004: Sujay Kumar; Initial Specification
!  03  Oct 2013: KR Arsenault;  Expanded to read in native STATSGO v1 files
!  10  Apr 2017: Shugong Wang; Read WHS soil parameter data and calculate
!                background soil albedo for JULES model 
! !INTERFACE:
subroutine create_soil_albedo_from_whs( n, nc_whs_param, soil_albedo) 

! !USES:
  !use LDT_paramDataMod ! check mask 
  use LDT_coreMod, only : LDT_rc
  use LDT_logMod,  only : LDT_logunit,  LDT_endrun
  use LDT_gridmappingMod    
  use LDT_paramDataMod
  use LDT_fileIOMod
  implicit none
! !ARGUMENTS: 
  integer,intent(in)    :: n
  character(len=128)    :: nc_whs_param    ! NetCDF file of CAP soil texture fraction file 
  real, intent(inout)   :: soil_albedo(LDT_rc%lnc(n),LDT_rc%lnr(n))
!
! !DESCRIPTION:
!  This subroutine read the soil texture fraction data and derive soil albedo parameters.
!  The fraction data is in lat-lon coordinate and at 1/20 degree resolution. 
!
!  The arguments are:
!  \begin{description}
!  \item[n]
!   index of the nest
!  \item[nc_whs_param]
!   file name of the NetCDF CAP soil texture fraction data 
!  \end{description}
!EOP
   
   integer, parameter :: input_cols = 7200
   integer, parameter :: input_rows = 3600
   real,    parameter :: input_xres = 1.0/20.0
   real,    parameter :: input_yres = 1.0/20.0

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
   real, allocatable :: subset_whs(:,:)  ! Read input parameter
   integer, allocatable  :: lat_line(:,:), lon_line(:,:)


   real, allocatable:: vegp_map(:,:), vegs_map(:,:), soil_code_map(:,:)
   real, allocatable:: whs_albedo(:,:) 
   real, allocatable:: modis_albedo(:,:) 
   real :: whs_table(34,2)
   integer :: soil_code
   real :: whs_code(23, 3) 
   
   whs_code = transpose(reshape((/0.0, 0.17, 0.25, &
                                 11.0, 0.26, 0.35, &
                                 12.0, 0.26, 0.35, &
                                 13.0, 0.26, 0.35, &
                                 14.0, 0.26, 0.35, &
                                 15.0, 0.26, 0.35, &
                                 16.0, 0.26, 0.35, &
                                 17.0, 0.17, 0.25, &
                                 18.0, 0.17, 0.25, &
                                 19.0, 0.17, 0.25, &
                                 20.0, 0.17, 0.25, &
                                 21.0, 0.17, 0.25, &
                                 22.0, 0.17, 0.25, &
                                 23.0, 0.11, 0.15, &
                                 24.0, 0.11, 0.15, &
                                 25.0, 0.11, 0.15, &
                                 26.0, 0.11, 0.15, &
                                 27.0, 0.11, 0.15, &
                                 28.0, 0.11, 0.15, &
                                 29.0, 0.26, 0.35, &
                                 30.0, 0.17, 0.25, &
                                 31.0, 0.11, 0.15, &
                                 34.0, 0.75, 0.75/), (/3,23/)))
   
   whs_table(1, 1) = whs_code(1,2)
   whs_table(1, 2) = whs_code(1,3)
   do i=2, 23
    soil_code = nint(whs_code(i, 1));
    whs_table(soil_code,1) = whs_code(i,2) 
    whs_table(soil_code,2) = whs_code(i,3) 
   enddo
! ___________________________________________________________________
   LDT_rc%waterclass = 0. 
   soil_albedo = LDT_rc%waterclass

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

   inquire(file=trim(nc_whs_param), exist=file_exists)
   if(.not.file_exists) then 
      write(LDT_logunit,*) "JULES WHS soil parameter file: ",trim(LDT_rc%txtfile(n))," not found."
      write(LDT_logunit,*) "Program stopping ..."
      call LDT_endrun
   endif

   write(unit=LDT_logunit,fmt=*) "[INFO] Reading JULES WHS soil parameter file: ",&
        trim(nc_whs_param)

   call LDT_checkDomainExtents(n, param_gridDesc(:))

   allocate(vegp_map(input_cols, input_rows)) 
   allocate(vegs_map(input_cols, input_rows)) 
   allocate(soil_code_map(input_cols, input_rows)) 
   allocate(whs_albedo(input_cols, input_rows)) 
   allocate(modis_albedo(input_cols, input_rows)) 

!- Read-in the JULSE CAP sand, silt and clay fraction data 
   call read_whs_param(nc_whs_param, vegp_map, vegs_map, soil_code_map, modis_albedo)

!   open(unit=1000, file='soil_code.txt', form='formatted', status='new', action='write' )
   whs_albedo(:,:) = LDT_rc%udef
   do c=1, input_cols
    do r=1, input_rows
      if(soil_code_map(c,r) > 0) then
        soil_code = nint(soil_code_map(c, r))
!        write(unit=1000, fmt='(I5)', advance='no') soil_code
        whs_albedo(c,r) = whs_table(soil_code, 1)
        if(vegp_map(c,r) .ge. 70 .and. vegp_map(c,r) .le. 73 .or. vegp_map(c,r) .eq. 36) then
          whs_albedo(c,r) = whs_table(soil_code, 2)
        endif
        if(vegs_map(c,r) .ge. 70 .and. vegs_map(c,r) .le. 73 .or. vegs_map(c,r) .eq. 36) then
          whs_albedo(c,r) = whs_table(soil_code, 2)
        endif
      endif
    enddo
   enddo

   ! replace MODIS missing data with WHS albedo
   do c=1, input_cols
    do r=1, input_rows
      if(modis_albedo(c,r) .eq. -1) then
        modis_albedo(c,r) = whs_albedo(c,r)
      endif
    enddo
  enddo
! -------------------------------------------------------------------
!    PREPARE SUBSETTED PARAMETER GRID FOR READING IN NEEDED DATA
! -------------------------------------------------------------------

!- Map Parameter Grid Info to LIS Target Grid/Projection Info -- 
   subparam_gridDesc = 0.
   call LDT_RunDomainPts( n, LDT_rc%soiltext_proj, param_gridDesc(:), &
            glpnc, glpnr, subpnc, subpnr, subparam_gridDesc, lat_line, lon_line )
   
   allocate( subset_whs(subpnc, subpnr) )
   subset_whs = LDT_rc%waterclass
!- Subset parameter read-in array:
   do r = 1, subpnr
      do c = 1, subpnc
         !subset_whs(c,r) = whs_albedo(lon_line(c,r), lat_line(c,r))
         subset_whs(c,r) = modis_albedo(lon_line(c,r), lat_line(c,r))
      enddo
   enddo

! -------------------------------------------------------------------
!     UPSCALING COARSE-SCALE GRIDS TO FINE LIS OUTPUT GRID
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
         gi(i) = subset_whs(c,r)
         if( gi(i) .ne. LDT_rc%udef ) li(i) = .true.
      enddo
   enddo
   ! need to do 


!- Transform parameter from original grid to LIS output grid:
   !LDT_rc%soiltext_gridtransform(n) = "mode"
   LDT_rc%soiltext_gridtransform(n) = "average"
   call LDT_transform_paramgrid(n, LDT_rc%soiltext_gridtransform(n), &
                      !param_gridDesc, mi, 1, gi, li, mo, go1, lo1 )
                      subparam_gridDesc, mi, 1, gi, li, mo, go1, lo1 )

!- Convert 1D count to 2D grid fgrds:
   i = 0
   do r = 1, LDT_rc%lnr(n)
      do c = 1, LDT_rc%lnc(n)
         i = i + 1
         soil_albedo(c,r) = go1(i)
        
         ! check mask 
         !if(LDT_LSMparam_struc(n)%landmask%value(c,r,1).eq.1 .and. soil_albedo(c,r) .eq. -9999) then
         ! LDT_LSMparam_struc(n)%landmask%value(c,r,1) = 0
         ! write(LDT_logunit,*) "Mask correction (soil albedo): col =  ", c, " row = ", r
         !endif
      enddo
   enddo
   

   deallocate( gi, li )


   deallocate(vegp_map) 
   deallocate(vegs_map) 
   deallocate(soil_code_map) 

   write(unit=LDT_logunit,fmt=*) "[INFO] Done creating JULES soil albedo based on WHS data"
 end subroutine create_soil_albedo_from_whs

 subroutine read_whs_param(nc_whs_param, vegp, vegs, soil_code, wsa)
   use netcdf
   implicit none 
   character(len=*) :: nc_whs_param
   real, intent(inout) :: vegp(7200,3600), soil_code(7200,3600), vegs(7200,3600), wsa(7200,3600)
   real, allocatable :: data_tmp(:,:) 
   integer :: ncid, varid
   
   ! open the file, NF90_NOWRITE tells netCDF we want read-only access to the file
   call jules_check_nc(nf90_open( trim(nc_whs_param), NF90_NOWRITE, ncid))
 
   call jules_check_nc(nf90_inq_varid(ncid, "whs_primary", varid))
   call jules_check_nc(nf90_get_var(ncid, varid, vegp))
   
   call jules_check_nc(nf90_inq_varid(ncid, "whs_secondary", varid))
   call jules_check_nc(nf90_get_var(ncid, varid, vegs))
 
   call jules_check_nc(nf90_inq_varid(ncid, "whs_soil_code", varid))
   call jules_check_nc(nf90_get_var(ncid, varid, soil_code))
   
   call jules_check_nc(nf90_inq_varid(ncid, "wsa_shortwave_soil_albedo", varid))
   call jules_check_nc(nf90_get_var(ncid, varid, wsa))
   
   call jules_check_nc(nf90_close(ncid))
 end subroutine read_whs_param
 
  
