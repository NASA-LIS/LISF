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
! !ROUTINE: Monfredaetal08_crops_module
!  \label{Monfredaetal08_crops_module}
!
! !REVISION HISTORY:
!  23  July 2014: KR Arsenault; Implemented new crop map features
!
! !INTERFACE:

module Monfredaetal08_crops_module
! !USES:
  use LDT_coreMod,  only : LDT_rc
  use LDT_logMod,   only : LDT_logunit, LDT_verify, LDT_endrun
  use LDT_gridmappingMod
  use LDT_fileIOMod
  use LDT_LSMCropModifier_Mod
#if ( defined USE_NETCDF3 || defined USE_NETCDF4 )
  use netcdf
#endif

 implicit none
 PRIVATE

!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  PUBLIC :: read_Monfredaetal08_croptype

CONTAINS

!BOP
!
! !ROUTINE: read_Monfredaetal08_croptype
!  \label{read_Monfredaetal08_croptype}
!
! !REVISION HISTORY:
!  23  July 2014: KR Arsenault; Implemented new crop map features
!
! !INTERFACE:
subroutine read_Monfredaetal08_croptype(n, num_types, fgrd)

! !USES:
  use LDT_coreMod,     only : LDT_rc
  use LDT_logMod,      only : LDT_logunit, LDT_getNextUnitNumber, &
       LDT_releaseUnitNumber, LDT_verify, LDT_endrun
  use LDT_gridmappingMod    
  use LDT_fileIOMod
  use LDT_paramTileInputMod, only: param_index_fgrdcalc
  use LDT_LSMCropModifier_Mod
#if ( defined USE_NETCDF3 || defined USE_NETCDF4 )
  use netcdf
#endif

  implicit none

! !ARGUMENTS: 
  integer, intent(in) :: n
  integer, intent(in) :: num_types
  real, intent(inout) :: fgrd(LDT_rc%lnc(n),LDT_rc%lnr(n),num_types)
!
! !DESCRIPTION:
!  This subroutine reads the Monfreda et al. (2008) crop data and returns the 
!  distribution of vegetation in each grid cell, in a lat/lon
!  projection.  
!
!  The arguments are:
!  \begin{description}
!   \item[n]
!     index of nest
!   \item[num_types]
!     number of landcover and/or just crop types
!   \item[fgrd]
!     fraction of grid covered by each vegetation and crop type
!   \end{description}
!EOP      

! Monfreda et al(2008) crop map:
   integer, parameter :: IN_cols = 4320
   integer, parameter :: IN_rows = 2160
   real,    parameter :: IN_xres = 1.0/12.0 ! or 360deg/4320 col points
   real,    parameter :: IN_yres = 1.0/12.0 ! or 180deg/2160 row points

! Other
   character(140) :: tempfile
   logical :: file_exists
   integer :: i, i2, j, t, c, r
   integer :: glpnc, glpnr             ! Parameter (global) total columns and rows
   integer :: subpnc, subpnr           ! Parameter subsetted columns and rows
   real    :: param_gridDesc(20)       ! Input parameter grid desc fgrd
   real    :: subparam_gridDesc(20)    ! Input parameter grid desc array
   integer, allocatable  :: lat_line(:,:), lon_line(:,:)

   character(20), allocatable :: croptype_array(:)
   real      :: cropdom(LDT_rc%lnc(n),LDT_rc%lnr(n))
   real      :: croptype_frac(LDT_rc%lnc(n),LDT_rc%lnr(n),LDT_rc%numcrop(n))
   real      :: subcroptype_frac(LDT_rc%lnc(n),LDT_rc%lnr(n),14)

!__________________________________________________________________

   croptype_frac = 0.0
   cropdom  = 0.0
   fgrd     = 0.0

!- Set parameter grid array inputs:
   LDT_LSMCrop_struc(n)%crop_proj = "latlon"
   param_gridDesc(1)  = 0.          ! Latlon
   param_gridDesc(2)  = IN_cols
   param_gridDesc(3)  = IN_rows
   param_gridDesc(4)  = -90.0  + (IN_yres/2) ! LL lat
   param_gridDesc(5)  = -180.0 + (IN_xres/2) ! LL lon
   param_gridDesc(6)  = 128
   param_gridDesc(7)  =  90.0 - (IN_yres/2)  ! UR lat
   param_gridDesc(8)  = 180.0 - (IN_xres/2)  ! UR lon
   param_gridDesc(9)  = IN_yres     ! dy: 0.083333
   param_gridDesc(10) = IN_xres     ! dx: 0.083333
   param_gridDesc(20) = 64

!- Map Parameter Grid Info to LIS Target Grid/Projection Info:
   subparam_gridDesc = 0.
   call LDT_RunDomainPts( n, LDT_LSMCrop_struc(n)%crop_proj, &
        param_gridDesc, glpnc, glpnr, subpnc, subpnr, &
        subparam_gridDesc, lat_line, lon_line )


!- Perform check on grid and projection choices selected:
   if( LDT_rc%lis_map_proj(n) == "latlon"   .or. &
       LDT_rc%lis_map_proj(n) == "mercator" .or. &
       LDT_rc%lis_map_proj(n) == "lambert" ) then
     if( param_gridDesc(10) .ne. (LDT_rc%gridDesc(n,9)/LDT_rc%lis_map_resfactor(n)) .and.&
         LDT_LSMCrop_struc(n)%crop_gridtransform .eq. "none" ) then
        write(LDT_logunit,*) "[ERR] 'Native' Monfreda et. al (2008) crop files have been selected, "
        write(LDT_logunit,*) "    with a resolution (0.0833deg), but the LIS run domain resolution"
        write(LDT_logunit,*) "    selected is not equal to that. So please select a spatial"
        write(LDT_logunit,*) "    transform other than 'none'."
        write(LDT_logunit,*) "Program stopping ..."
        call LDT_endrun
     endif
   endif

!- Check if Monfreda et al Crop map directory exists:
   tempfile = trim(LDT_LSMCrop_struc(n)%croptfile)//"/maize_5min.nc"
   inquire(file=tempfile, exist=file_exists)
   if(.not.file_exists) then
      write(LDT_logunit,*) "[ERR] The Monfreda et al. (2008) crop type map diretory: ",&
                            trim(LDT_LSMCrop_struc(n)%croptfile)," is not correct. "
      write(LDT_logunit,*) "Program stopping ..."
      call LDT_endrun
   endif
   write(LDT_logunit,*)"[INFO] Reading Monfreda et al. (2008) crop file: ",trim(LDT_LSMCrop_struc(n)%croptfile)
   write(LDT_logunit,*)" * Note: These crop files are continuous fields of crop areas, as in fractions." 

!- Retrieve the croptype array for user-selected classification:
   allocate( croptype_array(LDT_rc%numcrop(n)) )
   croptype_array = "none"
   call readcropinventory( n, LDT_LSMCrop_struc(n)%crop_classification, &
                           LDT_rc%numcrop(n), croptype_array ) 

!- Open each Crop-type file:
   do i = 1, LDT_rc%numcrop(n)

   !  Ignore "others" croptype of "CROPMAP" classification for now:
      if( croptype_array(i) == "others" ) cycle   

      tempfile = trim(LDT_LSMCrop_struc(n)%croptfile)//"/"//trim(croptype_array(i))//"_5min.nc"
      write(LDT_logunit,*)"  Opening crop type file: ",trim(tempfile)

      call readMonfredaCropfiles( n, tempfile, subpnc, subpnr, subparam_gridDesc, &
                                  lat_line, lon_line, croptype_frac(:,:,i) )

   end do    ! End Crop Type File Read


!- Some regions are dominated by "other" crop types, like California, Florida, 
!  if using the "CROPMAP" crop types....
   if( LDT_LSMCrop_struc(n)%crop_classification == "CROPMAP" ) then
     subcroptype_frac = 0.
     do i = 1, LDT_rc%numcrop(n)
        if( croptype_array(i) == "others" ) then
          do i2 = 1, 14
          ! California:
            if(i2==1) tempfile = trim(LDT_LSMCrop_struc(n)%croptfile)//"/rice_5min.nc"
            if(i2==2) tempfile = trim(LDT_LSMCrop_struc(n)%croptfile)//"/lettuce_5min.nc"
            if(i2==3) tempfile = trim(LDT_LSMCrop_struc(n)%croptfile)//"/grape_5min.nc"
            if(i2==4) tempfile = trim(LDT_LSMCrop_struc(n)%croptfile)//"/almond_5min.nc"
            if(i2==5) tempfile = trim(LDT_LSMCrop_struc(n)%croptfile)//"/maizefor_5min.nc"
            if(i2==6) tempfile = trim(LDT_LSMCrop_struc(n)%croptfile)//"/walnut_5min.nc"
            if(i2==7) tempfile = trim(LDT_LSMCrop_struc(n)%croptfile)//"/tomato_5min.nc"
            if(i2==8) tempfile = trim(LDT_LSMCrop_struc(n)%croptfile)//"/mixedgrass_5min.nc"
            if(i2==9) tempfile = trim(LDT_LSMCrop_struc(n)%croptfile)//"/alfalfa_5min.nc"
          ! Florida:
            if(i2==10) tempfile = trim(LDT_LSMCrop_struc(n)%croptfile)//"/orange_5min.nc"
            if(i2==11) tempfile = trim(LDT_LSMCrop_struc(n)%croptfile)//"/sugarcane_5min.nc"
            if(i2==12) tempfile = trim(LDT_LSMCrop_struc(n)%croptfile)//"/stringbean_5min.nc"
            if(i2==13) tempfile = trim(LDT_LSMCrop_struc(n)%croptfile)//"/greencorn_5min.nc"
            if(i2==14) tempfile = trim(LDT_LSMCrop_struc(n)%croptfile)//"/potato_5min.nc"
            write(LDT_logunit,*)"  Opening crop type file: ",trim(tempfile)

            call readMonfredaCropfiles( n, tempfile, subpnc, subpnr, subparam_gridDesc, &
                                        lat_line, lon_line, subcroptype_frac(:,:,i2) )
          enddo

       !- Find dominant crop type out of inventory list:
          do r = 1, LDT_rc%lnr(n)
             do c = 1, LDT_rc%lnc(n)
                if( maxval(subcroptype_frac(c,r,:)) == LDT_rc%udef ) then
                   croptype_frac(c,r,i) = LDT_rc%udef
                else
                   croptype_frac(c,r,i) = maxval(subcroptype_frac(c,r,:),1,&
                                        mask=subcroptype_frac(c,r,:).ne.LDT_rc%udef )
                endif
             enddo
          enddo
          exit   ! Exit main croptype loop when "others" reached"
        endif    
     enddo
   endif         ! End accounting for CROPMAP "others" crop type


!- Find dominant crop type out of inventory list:
   do r = 1, LDT_rc%lnr(n)
      do c = 1, LDT_rc%lnc(n)
         cropdom(c,r) = maxval(croptype_frac(c,r,:) )

         if( maxval(croptype_frac(c,r,:)) == LDT_rc%udef ) then
            fgrd(c,r,1) = LDT_rc%udef
         else
            if( LDT_LSMCrop_struc(n)%crop_classification == "CROPMAP" ) then
               fgrd(c,r,1) = maxloc(croptype_frac(c,r,:),1,&
                             mask=croptype_frac(c,r,:).ne.LDT_rc%udef )+13.0
!               if(fgrd(c,r,1) == 21.) fgrd(c,r,1) = LDT_rc%udef
            else
              fgrd(c,r,1) = maxloc(croptype_frac(c,r,:),1,&
                            mask=croptype_frac(c,r,:).ne.LDT_rc%udef )
            endif
         endif
      enddo
   enddo

   deallocate( croptype_array )
   deallocate( lat_line, lon_line )

  write(LDT_logunit,*) "[INFO] Done reading 'Monfredaetal08' Crop map type. "


end subroutine read_Monfredaetal08_croptype


subroutine readMonfredaCropfiles( n, tempfile, subpnc, subpnr, subparam_gridDesc, &
                                  lat_line, lon_line, croptype_frac )

   integer,       intent(in) :: n
   character(140),intent(in) :: tempfile
   integer,       intent(in) :: subpnc, subpnr    
   real,          intent(in) :: subparam_gridDesc(20)  
   integer,       intent(in) :: lat_line(subpnc,subpnr), lon_line(subpnc,subpnr)
   real,         intent(out) :: croptype_frac(LDT_rc%lnc(n),LDT_rc%lnr(n))

! Netcdf file read-in entries:
   integer :: cropid
   integer :: latid, lonid
   integer :: nrows, ncols
! Other
   logical :: file_exists
   integer :: ftn, ierr
   integer :: j, t, c, r
   integer :: mi                       ! Total number of input param grid array points
   integer :: mo                       ! Total number of output LIS grid array points
   real,    allocatable  :: gi1(:)     ! Input parameter 1d grid
   logical*1,allocatable :: li1(:)     ! Input logical mask (to match gi)
   real      :: go1(LDT_rc%lnc(n)*LDT_rc%lnr(n))            ! Output lis 1d grid
   logical*1 :: lo1(LDT_rc%lnc(n)*LDT_rc%lnr(n))            ! Output logical mask (to match go)

   real, allocatable :: read_file(:,:)        ! Read in parameter
   real, allocatable :: yrev_file(:,:)        ! Y-reverse parameter
   real, allocatable :: read_crop_sub(:,:)    ! Read input parameter

!_______________________________________________________________________________________

   if( allocated(read_file) ) then
       deallocate(read_file)
   endif
!
!  READ IN MONFREDA ETAL (2008) CROP TYPE PARAMETER FIELDS 
!
#if ( defined USE_NETCDF3 || defined USE_NETCDF4 )
    ierr = nf90_open(path=tempfile,mode=NF90_NOWRITE,ncid=ftn)
    call LDT_verify(ierr,'error opening Monfreda etal (2008) crop type file')

    ierr = nf90_inq_varid(ftn,'cropdata',cropid)
    call LDT_verify(ierr, 'nf90_inq_varid failed for crop type in read_Monfredaetal08_croptype')

    ierr = nf90_inq_dimid(ftn,'longitude',lonid)
    call LDT_verify(ierr,'nf90_inq_dimid failed for longitude in read_Monfredaetal08_croptype')
    ierr = nf90_inq_dimid(ftn,'latitude',latid)
    call LDT_verify(ierr,'nf90_inq_dimid failed for latitude in read_Monfredaetal08_croptype')
    ierr = nf90_inquire_dimension(ftn,lonid,len=ncols)
    call LDT_verify(ierr,'nf90_inquire_dimension for longitude')
    ierr = nf90_inquire_dimension(ftn,latid,len=nrows)
    call LDT_verify(ierr,'nf90_inquire_dimension for latitude')

    allocate( read_file(ncols,nrows) )
    ierr = nf90_get_var(ftn, cropid, read_file)
    call LDT_verify(ierr, 'nf90_get_var failed for Monfreda etal(2008) crop type')

    ierr = nf90_close(ftn)
    call LDT_verify(ierr, 'nf90_close failed in read_Monfredaetal08_croptype')
#endif
 !- Reverse-Y and Convert 8-bit unsigned integers:
    allocate( yrev_file(ncols,nrows) )
    yrev_file = LDT_rc%udef
    j = 0
    do r = nrows, 1, -1
       j = j + 1
       do c = 1, ncols
          yrev_file(c,j) = read_file(c,r)
       end do
    end do
    deallocate( read_file )
 !- Subset parameter read-in array:
    allocate( read_crop_sub(subpnc, subpnr) )
    read_crop_sub = LDT_rc%udef
    do r = 1, subpnr
       do c = 1, subpnc
          read_crop_sub(c,r) = yrev_file(lon_line(c,r),lat_line(c,r))
       enddo
    enddo
    deallocate(yrev_file)

! -------------------------------------------------------------------
!   UPSCALING/DOWNSCALING GRIDS TO LIS OUTPUT GRID
! -------------------------------------------------------------------
    mi = nrows*ncols
    allocate( gi1(mi), li1(mi) )
    gi1 = LDT_rc%udef;  li1 = .false.
    mo = LDT_rc%lnc(n)*LDT_rc%lnr(n)
    lo1 = .false.;  go1 = 0.
 !- Assign 2-D array to 1-D for aggregation routines:
    j = 0
    do r = 1, subpnr
       do c = 1, subpnc;  j = j + 1
          gi1(j) = read_crop_sub(c,r)
          if( gi1(j) <= 0. .or. gi1(j) > 100.0 ) gi1(j) = LDT_rc%udef
          if( gi1(j) .ne. LDT_rc%udef )  li1(j) = .true.
       enddo
    enddo
    deallocate( read_crop_sub )

 !- Upscale/downscale routines:
    select case ( LDT_LSMCrop_struc(n)%crop_gridtransform )

   !- Single layer transformation:
      case( "none", "average", "neighbor", "bilinear", "budget-bilinear" )

       !- Transform parameter from original grid to LIS output grid:
          call LDT_transform_paramgrid(n, LDT_LSMCrop_struc(n)%crop_gridtransform, &
                   subparam_gridDesc, mi, 1, gi1, li1, mo, go1, lo1 )

       !- Convert 1D soil fractions to 2D grid arrays:
          j = 0
          do r = 1, LDT_rc%lnr(n)
             do c = 1, LDT_rc%lnc(n)
                j = j + 1
                croptype_frac(c,r) = go1(j)
              ! Make sure all negative values are set to universal undefined value:
                if( croptype_frac(c,r) < 0.  ) croptype_frac(c,r) = LDT_rc%udef
                if( croptype_frac(c,r) > 100.) croptype_frac(c,r) = LDT_rc%udef
             enddo
          enddo

     case default
       write(*,*)"[WARN] Other spatial grid transformations are not currently supported"
       write(*,*)"   for the tiled 'Monfredaetal08' crop type maps. Please select either:"
       write(*,*)"  -- neighbor, bilinear, budget-bilinear (to downscale)"
       write(*,*)"  -- average (to upscale)"
       write(*,*)" Stopping program ..."
       call LDT_endrun

    end select  ! End vegtype/cnt aggregation method
    deallocate( gi1, li1 )

 end subroutine readMonfredaCropfiles

end module Monfredaetal08_crops_module
