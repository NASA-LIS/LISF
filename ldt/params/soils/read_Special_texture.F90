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
! !ROUTINE: read_Special_texture
! \label{read_Special_texture}
!
! !REVISION HISTORY:
!  03 Nov 2012: KR Arsenault; Initial Specification
!  03 May 2014: KR Arsenault; Updated for global files
!
! !INTERFACE:
subroutine read_Special_texture( n, num_bins, array, texture_layers )

! !USES:
  use LDT_coreMod,   only : LDT_rc
  use LDT_logMod,    only : LDT_logunit, LDT_verify, LDT_endrun
  use LDT_fileIOMod, only : LDT_checkDomainExtents, LDT_transform_paramgrid
  use LDT_paramTileInputMod, only: param_index_fgrdcalc
  use LDT_gridmappingMod
#if ( defined USE_NETCDF3 || defined USE_NETCDF4 )
  use netcdf
#endif

  implicit none
! !ARGUMENTS: 
  integer,intent(in)    :: n
  integer,intent(in)    :: num_bins   ! Number of soil types
  real,   intent(inout) :: array(LDT_rc%lnc(n),LDT_rc%lnr(n),num_bins)
  real,   intent(inout) :: texture_layers(LDT_rc%lnc(n),LDT_rc%lnr(n),1)
!
! !DESCRIPTION:
!  This subroutine retrieves Special soil texture data and reprojects
!  it to the latlon projection. 
!
!  The arguments are:
!  \begin{description}
!  \item[n]
!   index of the nest
!  \item[num_bins]
!   number of soil types
!  \item[array]
!   output field with the retrieved soil texture
!  \end{description}
!EOP
  integer :: c, r, t, i
  integer :: ftn, ierr
  logical :: file_exists
  integer :: water_class
! Netcdf file read-in entries:
  integer :: textureid
  integer :: latid, lonid
  integer :: nrows, ncols

!- Used for upscaling:
  integer :: glpnc, glpnr             ! Parameter (global) total columns and rows
  integer :: subpnc, subpnr           ! Parameter subsetted columns and rows
  integer :: mi                       ! Total number of input param grid array points
  integer :: mo                       ! Total number of output LIS grid array points
  real    :: param_gridDesc(20)       ! Input parameter grid desc array
  real    :: subparam_gridDesc(20)    ! Subsetted Input parameter grid desc array
  integer, allocatable :: lat_line(:,:), lon_line(:,:)
!  integer, allocatable :: n11(:)                ! Maps each input grid point to output grid.
  real,    allocatable :: gi1(:)                ! input parameter 1d grid
  logical*1,allocatable:: li1(:)                ! input logical mask (to match gi)
  real      :: go1(LDT_rc%lnc(n)*LDT_rc%lnr(n)) ! output lis 1d grid
  logical*1 :: lo1(LDT_rc%lnc(n)*LDT_rc%lnr(n)) ! output logical mask (to match go)
  real      :: go2(LDT_rc%lnc(n)*LDT_rc%lnr(n), LDT_rc%nt) ! Output lis 1d grid
  logical*1 :: lo2(LDT_rc%lnc(n)*LDT_rc%lnr(n), LDT_rc%nt) ! Output logical mask (to match go)

  real    :: soiltext(LDT_rc%lnc(n),LDT_rc%lnr(n))
  real    :: soilcnt(LDT_rc%lnc(n),LDT_rc%lnr(n),num_bins)
  real    :: fgrd(LDT_rc%lnc(n),LDT_rc%lnr(n),num_bins)
  real, allocatable :: read_texture(:,:)        ! Read in parameter
  real, allocatable :: yrev_texture(:,:)        ! Y-reverse parameter
  real, allocatable :: subset_texture(:,:)      ! Subset input parameter

! ___________________________________________________________________

  water_class = 14    ! Water class for USDA types
  
  soiltext = LDT_rc%udef
  fgrd = 0.
  array = 0.
  texture_layers = 0.

  write(unit=LDT_logunit,fmt=*) "MSG: Reading Special texture file: ",&
                                 trim(LDT_rc%txtfile(n))

  inquire(file=trim(LDT_rc%txtfile(n)), exist=file_exists)
  if(.not.file_exists) then 
     write(LDT_logunit,*) "Texture map ",trim(LDT_rc%txtfile(n))," not found."
     write(LDT_logunit,*) "Program stopping ..."
     call LDT_endrun
  endif

#if ( defined USE_NETCDF3 || defined USE_NETCDF4 )

  ierr = nf90_open(path=trim(LDT_rc%txtfile(n)),mode=NF90_NOWRITE,ncid=ftn)
  call LDT_verify(ierr,'error opening USDA soil texture data')

  ierr = nf90_inq_varid(ftn,'Soil_type',textureid)
!  ierr = nf90_inq_varid(ftn,'sand_percent',textureid)
  call LDT_verify(ierr, 'nf90_inq_varid failed for texture in read_Special_texture')

  ierr = nf90_inq_dimid(ftn,'lon',lonid)
  call LDT_verify(ierr,'nf90_inq_dimid failed for longitude in read_Special_texture')

  ierr = nf90_inq_dimid(ftn,'lat',latid)
  call LDT_verify(ierr,'nf90_inq_dimid failed for latitude in read_Special_texture')

  ierr = nf90_inquire_dimension(ftn,lonid,len=ncols)
  call LDT_verify(ierr,'nf90_inquire_dimension for longitude')

  ierr = nf90_inquire_dimension(ftn,latid,len=nrows)
  call LDT_verify(ierr,'nf90_inquire_dimension for latitude')

  allocate( read_texture(ncols,nrows) )
  ierr = nf90_get_var(ftn, textureid, read_texture)
  call LDT_verify(ierr, 'nf90_get_var failed for texture')

  ierr = nf90_close(ftn)
  call LDT_verify(ierr, 'nf90_close failed in read_Special_texture')

#endif
!____________________________________________________________________________________

 !- Set parameter grid array inputs:
    param_gridDesc(1)  = 0.    ! Latlon
    param_gridDesc(2)  = float(ncols)
    param_gridDesc(3)  = float(nrows)
    param_gridDesc(4)  = -60.0000   ! LL lat
    param_gridDesc(5)  = -180.0000  ! LL lon
    param_gridDesc(6)  = 128
    param_gridDesc(7)  = 90.0000    ! UR lat
    param_gridDesc(8)  = 180.0000   ! UR lon
    param_gridDesc(9)  = 0.009000900090009
    param_gridDesc(10) = 0.009000900090009
    param_gridDesc(20) = 64

! -------------------------------------------------------------------
!    PREPARE SUBSETTED PARAMETER GRID FOR READING IN NEEDED DATA
! -------------------------------------------------------------------

 ! Reverse-Y and Convert 8-bit unsigned integers:
   allocate( yrev_texture(ncols,nrows) )
   yrev_texture = water_class
   i = 0
   do r = nrows, 1, -1
      i = i + 1
      do c = 1, ncols
         yrev_texture(c,i) = read_texture(c,r)
      end do
   end do
   deallocate( read_texture )

!- Map Parameter Grid Info to LIS Target Grid/Projection Info -- 
   subparam_gridDesc = 0.
   call LDT_RunDomainPts( n, LDT_rc%lc_proj, param_gridDesc(:), &
            glpnc, glpnr, subpnc, subpnr, subparam_gridDesc, lat_line, lon_line )

   allocate( subset_texture(subpnc, subpnr) )
   subset_texture = water_class

!- Subset parameter read-in array:
   do r = 1, subpnr
      do c = 1, subpnc
         subset_texture(c,r) = yrev_texture(lon_line(c,r),lat_line(c,r))
      enddo
   enddo
   deallocate(yrev_texture)


! -------------------------------------------------------------------
!    UPSCALING/DOWNSCALING GRIDS TO LIS OUTPUT GRID
! -------------------------------------------------------------------
   mi = nrows*ncols
   allocate( gi1(mi), li1(mi) )! , n11(mi) )
   gi1 = water_class
   li1 = .false.

   mo = LDT_rc%lnc(n)*LDT_rc%lnr(n)
   lo1 = .false.

!- Assign 2-D array to 1-D for aggregation routines:
   i = 0
   do r = 1, subpnr
      do c = 1, subpnc;  i = i + 1
         gi1(i) = subset_texture(c,r)
         if( gi1(i) == 0.    )  gi1(i) = 13.                  ! Set gravel class (=0) to 13
         if( gi1(i) == -128. )  gi1(i) = float(water_class)   ! Set water (=-128) to 14
         if( gi1(i) ==  127. )  gi1(i) = float(water_class)   ! Set water (=-128) to 14
         if( gi1(i) .ne. LDT_rc%udef )  li1(i) = .true.
      enddo
   enddo
   deallocate( subset_texture )

!- Aggregation/Spatial Transform Section:
   select case ( LDT_rc%soiltext_gridtransform(n) )

  !- (a) Estimate NON-TILED dominant soil texture types:
     case( "neighbor", "mode" )

    !- Transform parameter from original grid to LIS output grid:
       call LDT_transform_paramgrid(n, LDT_rc%soiltext_gridtransform(n), &
                subparam_gridDesc, mi, 1, gi1, li1, mo, go1, lo1 )

    !- Convert 1D soiltext to 2D grid arrays:
       i = 0
       do r = 1, LDT_rc%lnr(n)
          do c = 1, LDT_rc%lnc(n)
             i = i + 1
             soiltext(c,r) = go1(i)
          enddo
       enddo

  !- (b) Estimate TILED land cover files (soiltext):
     case( "tile" )

     !- Transform parameter from original grid to LIS output grid:
        call LDT_transform_paramgrid(n, LDT_rc%soiltext_gridtransform(n), &
                 subparam_gridDesc, mi, LDT_rc%nt, gi1, li1, mo, go2, lo2 )

     !- Convert 1D soiltext to 2D grid arrays:
        i = 0
        do r = 1, LDT_rc%lnr(n)
           do c = 1, LDT_rc%lnc(n)
              i = i + 1
              do t = 1, LDT_rc%nt
                 soilcnt(c,r,t) = go2(i,t)
              end do
           enddo
        enddo
        array = soilcnt 

     case default
       write(LDT_logunit,*) "Soil texture spatial transform option,",&
             trim(LDT_rc%soiltext_gridtransform(n))," not available."
       write(LDT_logunit,*) "Program stopping ..."
       call LDT_endrun

   end select  ! End soiltype/cnt aggregation method
   deallocate( gi1, li1 )


!- Bring 2-D Array to 3-D Soil tile space:
   if( LDT_rc%soiltext_gridtransform(n) == "none"  .or. &
       LDT_rc%soiltext_gridtransform(n) == "mode"  .or. &
       LDT_rc%soiltext_gridtransform(n) == "neighbor" ) then

      do r = 1, LDT_rc%lnr(n)
         do c = 1, LDT_rc%lnc(n)
            if( nint(soiltext(c,r)) .le. 0 ) then
               soiltext(c,r) = float(water_class)
            endif
            if ((nint(soiltext(c,r)) .ne. water_class ) .and. &
                (nint(soiltext(c,r)) .ne. LDT_rc%udef)) then
               array(c,r,NINT(soiltext(c,r))) = 1.0
          endif
         enddo
      enddo
   end if

!- Estimate fraction of grid (fgrid) represented by soil type::
   call param_index_fgrdcalc( n, LDT_rc%soiltext_proj, LDT_rc%soiltext_gridtransform(n), &
                              water_class, num_bins, array, fgrd )
   array = fgrd

   write(unit=LDT_logunit,fmt=*) "MSG: Done reading Special texture file."

 end subroutine read_Special_texture

