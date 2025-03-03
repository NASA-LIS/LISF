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
! !ROUTINE: read_um_soil_params
!  \label{read_um_soil_params}
!
! !REVISION HISTORY:
!  03 Sept 2004: Sujay Kumar; Initial Specification
!  23  Aug 2013: KR Arsenault; Implemented new features to read different
!                 resolutions and generate landmask from landcover
!  23  Apr 2014: KR Arsenault;  Added new optimized interpolation code
!  23  Mar 2017: Shugong Wang; Modified for creating PFT for JULES Modified IGBP
!  22  Jun 2017: Shugong Wang; Modified for creating PFT for JULES based on UKMO IGBP
!  11  Sep 2017: Shugong Wang; Modified for creating PFT for JULES based on UM Ancillary 
!  13  Sep 2017: Shugong Wang; Modified for reading UM soil aprameters 
! !INTERFACE:
subroutine read_um_soil_params(n, nc_file, nc_var, param_array)

! !USES:
  use LDT_coreMod,     only : LDT_rc
  use LDT_logMod,      only : LDT_logunit, LDT_getNextUnitNumber, &
             LDT_releaseUnitNumber, LDT_verify, LDT_endrun
  use LDT_gridmappingMod    
  use LDT_fileIOMod
  use LDT_paramTileInputMod, only: param_index_fgrdcalc_pft
  use netcdf 
  implicit none

! !ARGUMENTS: 
  integer, intent(in) :: n
  character(len=*), intent(in)   :: nc_file
  character(len=*), intent(in)   :: nc_var 
  real, intent(inout) :: param_array(LDT_rc%lnc(n),LDT_rc%lnr(n))

!
! !DESCRIPTION:
!  This subroutine reads the MODIS landcover data and returns the 
!  distribution of vegetation in each grid cell, in a lat/lon
!  projection.  Also, the landmask is either generated and/or 
!  read in this routine.
!
!  The arguments are:
!  \begin{description}
!   \item[n]
!     index of nest
!   \item[fgrd]
!     fraction of grid covered by each vegetation type
!   \item[param_array]
!     soil parameters for the region of interest
!   \end{description}
!EOP      
!
! IGBP-NCEP landcover version:
    integer, parameter :: IN_cols_10k = 2560  !43200 
    integer, parameter :: IN_rows_10k = 1920  !21600
    integer, parameter :: IN_cols_1k = 25600  !43200 
    integer, parameter :: IN_rows_1k = 19200  !21600
    real,    parameter :: IN_xres = 360.0/25600.0
    real,    parameter :: IN_yres = 180.0/19200.0 
    real, allocatable :: param_1k(:,:)
    real, allocatable :: param_10k(:,:)

    integer :: ftn, ierr, ios1
    logical :: file_exists
    integer :: i, t, c, r, r10k, c10k, npft
    integer :: input_cols, input_rows
    integer :: glpnc, glpnr                                  ! Parameter (global) total columns and rows
    integer :: subpnc, subpnr                                ! Parameter subsetted columns and rows
    integer :: mi                                            ! Total number of input param grid array points
    integer :: mo                                            ! Total number of output LIS grid array points
    real    :: param_gridDesc(20)                            ! Input parameter grid desc array
    real    :: subparam_gridDesc(20)                         ! Subsetted parameter grid desc array
    integer, allocatable  :: lat_line(:,:), lon_line(:,:)
    real,    allocatable  :: gi(:)                           ! Input parameter 1d grid
    logical*1,allocatable :: li(:)                           ! Input logical mask (to match gi)
    real      :: go1(LDT_rc%lnc(n)*LDT_rc%lnr(n))            ! Output lis 1d grid
    logical*1 :: lo1(LDT_rc%lnc(n)*LDT_rc%lnr(n))            ! Output logical mask (to match go)

    real, allocatable :: subset_param(:,:)                   ! Read input parameter
    integer   :: ncid, varid, ncstatus 
    character*50 :: grid_transform_mode

    allocate(param_1k(IN_cols_1k,IN_rows_1k))
    allocate(param_10k(IN_cols_10k,IN_rows_10k))


!- Check if land cover file exists:
    !inquire( file=trim(LDT_rc%vfile(n)), exist=file_exists )
    inquire( file=trim(nc_file), exist=file_exists )    
    if(.not. file_exists) then
       write(LDT_logunit,*) "The UM soil parameter NetCDF file: ",trim(nc_file)," does not exist. "
       write(LDT_logunit,*) "Program stopping ..."
       call LDT_endrun
    endif

    ftn = LDT_getNextUnitNumber()
    !write(LDT_logunit,*) "[INFO] Reading UM soil parameter file: ",trim(LDT_rc%vfile(n))
    write(LDT_logunit,*) "[INFO] Reading UM soil parameter file: ",trim(nc_file)

    input_cols = IN_cols_1k
    input_rows = IN_rows_1k

!- Set parameter grid array inputs:
    param_gridDesc(1)  = 0.          ! Latlon
    param_gridDesc(2)  = input_cols
    param_gridDesc(3)  = input_rows
    param_gridDesc(4)  = -90.0  + (IN_yres/2) ! LL lat
    param_gridDesc(5)  = -180.0 + (IN_xres/2) ! LL lon
    param_gridDesc(6)  = 128
    param_gridDesc(7)  =  90.0 - (IN_yres/2)  ! UR lat
    param_gridDesc(8)  = 180.0 - (IN_xres/2)  ! UR lon
    param_gridDesc(9)  = IN_xres     ! 
    param_gridDesc(10) = IN_yres     !
    param_gridDesc(20) = 64

! ! Open file:

!    ncstatus = nf90_open(LDT_rc%vfile(n), NF90_NOWRITE, ncid)
    ncstatus = nf90_open(nc_file, NF90_NOWRITE, ncid)    
    if (ncstatus /= nf90_noerr) write(LDT_logunit,*) "[INFO] cannot open ", trim(nc_file)
    ncstatus = nf90_inq_varid(ncid, trim(nc_var), varid)
    if (ncstatus /= nf90_noerr) write(LDT_logunit,*) "[INFO] UM soil parameter ", trim(nc_var), " not exist in ",  trim(nc_file)
    ncstatus = nf90_get_var(ncid, varid, param_10k)
    if (ncstatus /= nf90_noerr) write(LDT_logunit,*) "[INFO] cannot access dataset ", trim(nc_var), " in ",  trim(nc_file)
  
! ! downscale 10km PFT into 1km PFT
    do c=1, IN_cols_1k
      c10k = ceiling(c/10.0)
      do r=1, IN_rows_1k
        r10k = ceiling(r/10.0)
        if((param_10k(c10k,r10k) .ge. -1e10) .and. (param_10k(c10k,r10k) .le. 1e10)) then
          param_1k(c, r) = param_10k(c10k,r10k)
        else
          param_1k(c, r) = LDT_rc%udef
        endif
      enddo
    enddo
    
    write(LDT_logunit,*) "[INFO] Done reading ", trim(nc_file) , " for data set ", trim(nc_var) 


! -------------------------------------------------------------------
!    PREPARE SUBSETTED PARAMETER GRID FOR READING IN NEEDED DATA
!  - Map Parameter Grid Info to LIS Target Grid/Projection Info -- 
    subparam_gridDesc = 0.
    call LDT_RunDomainPts( n, LDT_rc%lc_proj, param_gridDesc(:), &
             glpnc, glpnr, subpnc, subpnr, subparam_gridDesc, lat_line, lon_line )

    allocate( subset_param(subpnc, subpnr) )
    subset_param = LDT_rc%udef

!- Subset parameter read-in array:
    do r = 1, subpnr
       do c = 1, subpnc
          subset_param(c, r) = param_1k(lon_line(c,r),lat_line(c, r))
       enddo
    enddo

! -------------------------------------------------------------------
!     AGGREGATING FINE-SCALE GRIDS TO COARSER LIS OUTPUT GRID
! -------------------------------------------------------------------
    mi = subpnc*subpnr
    mo = LDT_rc%lnc(n)*LDT_rc%lnr(n)
    allocate( gi(mi), li(mi) )
    
    gi = LDT_rc%udef
    li = .false.
    lo1 = .false.

    !- Assign 2-D array to 1-D for aggregation routines:
    i = 0
    do r = 1, subpnr
       do c = 1, subpnc
          i = i + 1
          gi(i) = subset_param(c,r)
          if( gi(i) .ne. LDT_rc%udef ) li(i) = .true.
       enddo
    enddo
    
    grid_transform_mode="neighbor"
    !- call LDT_transform_paramgrid(n, LDT_rc%lc_gridtransform(n), &
    call LDT_transform_paramgrid(n, grid_transform_mode, &
               subparam_gridDesc, mi, 1, gi, li, mo, go1, lo1 )

    !- Convert 1D vegcnt to 2D grid arrays:
    i = 0
    do r = 1, LDT_rc%lnr(n)
      do c = 1, LDT_rc%lnc(n)
        i = i + 1
        param_array(c, r) = go1(i)
      enddo
    enddo


    ! free memory 
    deallocate( gi, li )
    deallocate( subset_param )
    deallocate( param_1k)
    deallocate( param_10k)

    call LDT_releaseUnitNumber(ftn)

end subroutine read_um_soil_params
