!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA GSFC Land Data Toolkit (LDT) V1.0
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LDT_misc.h"
!BOP
!
! !ROUTINE: MIRCA2000_crops_module
!  \label{MIRCA2000_crops_module}
!
! !REVISION HISTORY:
!  21  May 2019: H Beaudoing; Implemented MIRCA2000 crop map features
!                             The routine is based on ReadProcess_MIRCA 
!                             section of preprocess_irrig.F90, by Sarith
!
! !INTERFACE:

module MIRCA2000_crops_module
! !USES:
  use LDT_coreMod,  only : LDT_rc
  use LDT_logMod,   only : LDT_logunit, LDT_verify, LDT_endrun
  use LDT_gridmappingMod
  use LDT_fileIOMod
  use LDT_LSMCropModifier_Mod

 implicit none
 PRIVATE

!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  PUBLIC :: read_MIRCA2000_croptype

CONTAINS

!BOP
!
! !ROUTINE: read_MIRCA2000_croptype
!  \label{read_MIRCA2000_croptype}
!
! !REVISION HISTORY:
!  21  May 2019: H Beaudoing; Implemented MIRCA2000 crop map features
!
! !INTERFACE:
subroutine read_MIRCA2000_croptype(n, num_types, fgrd)

! !USES:
  use LDT_coreMod,     only : LDT_rc
  use LDT_logMod,      only : LDT_logunit, LDT_getNextUnitNumber, &
       LDT_releaseUnitNumber, LDT_verify, LDT_endrun
  use LDT_gridmappingMod    
  use LDT_fileIOMod
  use LDT_paramTileInputMod, only: param_index_fgrdcalc
  use LDT_LSMCropModifier_Mod

  implicit none

! !ARGUMENTS: 
  integer, intent(in) :: n
  integer, intent(in) :: num_types
  real, intent(inout) :: fgrd(LDT_rc%lnc(n),LDT_rc%lnr(n),num_types)
!
! !DESCRIPTION:
!  This subroutine reads the MIRCA2000 Portmann et al. (2010) crop data 
!  and returns the distribution of vegetation in each grid cell, in a lat/lon
!  projection.  
!  The crop fraction inputs are avaiable at monthly. Currently, the monthly
!  information is condensed to annual values in the reading routine.
!
!  Number of crops for MIRCA2000 can be specified in ldt.config 
!  (via LDT_rc%numcrop) to either 26 or 52 where:
!  LDT_rc%numcrop = 26 => irrigated crops only
!  LDT_rc%numcrop = 52 => irrigated + rainfed crops 
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

   integer, parameter :: IN_NCROPS = 26 ! crop types, repeated for irrigated and                                        ! rainfed

! Other
   character(140) :: tempfile
   logical :: file_exists
   integer :: i, ii, j, t, c, r

   character(20), allocatable :: croptype_array(:)
   real      :: cropdom(LDT_rc%lnc(n),LDT_rc%lnr(n))
   real      :: croptype_frac(LDT_rc%lnc(n),LDT_rc%lnr(n),LDT_rc%numcrop(n))
   real   :: crop_array(LDT_rc%numcrop(n))

!__________________________________________________________________

   croptype_frac = 0.0
   cropdom  = 0.0
   fgrd     = 0.0
   crop_array = 0.0

!- Check if MIRCA2000 Crop map directory exists:
   tempfile = trim(LDT_LSMCrop_struc(n)%croptfile)//"01_irrigated_12.flt"
   inquire(file=tempfile, exist=file_exists)
   if(.not.file_exists) then
      write(LDT_logunit,*) "[ERR] The MIRCA2000 crop type map diretory: ",&
                            trim(LDT_LSMCrop_struc(n)%croptfile)," is not correct. "
      write(LDT_logunit,*) "Program stopping ..."
      call LDT_endrun
   endif
   write(LDT_logunit,*)"[INFO] Reading MIRCA2000 crop file: ",trim(LDT_LSMCrop_struc(n)%croptfile)
   write(LDT_logunit,*)" * Note: These crop files are continuous fields of crop areas (ha), as in fractions." 

!- Retrieve the croptype array for user-selected classification:
   allocate( croptype_array(LDT_rc%numcrop(n)) )
   croptype_array = "none"
   call readcropinventory( n, LDT_rc%crop_classification(n), &
                           LDT_rc%numcrop(n), croptype_array ) 

!- Open each Crop-type file:
!croptfile = '/gpfsm/dnb31/hkato/IRRIGATION/MIRCA2000/monthly_growing_area_grids/crop_'

   do i = 1, LDT_rc%numcrop(n)

    if ( i .le. IN_NCROPS ) then   ! i=1-26
      tempfile = trim(LDT_LSMCrop_struc(n)%croptfile)//trim(croptype_array(i))//"_irrigated_12.flt"
    else   ! i= 27-52
      ii = LDT_rc%numcrop(n) - i
      tempfile = trim(LDT_LSMCrop_struc(n)%croptfile)//trim(croptype_array(ii))//"_rainfed_12.flt"
    endif
    write(LDT_logunit,*)"  Opening crop type file: ",trim(tempfile)

    call readMIRCA2000Cropfiles( n, tempfile, croptype_frac(:,:,i) )
                                 

   end do    ! End Crop Type File Read

!- Write final gridcell fractions:
   do r = 1, LDT_rc%lnr(n)
      do c = 1, LDT_rc%lnc(n)
!- Find dominant crop type out of inventory list:
         crop_array(:) = croptype_frac(c,r,:)
         cropdom(c,r) = maxloc(crop_array,1)

         if( maxval(croptype_frac(c,r,:)) == LDT_rc%udef ) then
            fgrd(c,r,:) = LDT_rc%udef
         else
           do i = 1, LDT_rc%numcrop(n)
            fgrd(c,r,i) = croptype_frac(c,r,i)
           enddo
         endif
      enddo
   enddo

   deallocate( croptype_array )

  write(LDT_logunit,*) "[INFO] Done reading 'MIRCA2000' Crop map type. "

!HKB compute irrigfrac, rainfedfrac, and paddy??

end subroutine read_MIRCA2000_croptype


 subroutine readMIRCA2000Cropfiles( n, tempfile, croptype_frac )
                                  
! !USES:
    use LDT_coreMod
    use LDT_gridmappingMod
    use LDT_paramTileInputMod
    use LDT_LSMCropModifier_Mod
    use LDT_fileIOMod
    use LDT_logMod

    implicit none
! !ARGUMENTS:
   integer,       intent(in) :: n
   character(140),intent(in) :: tempfile
   real,         intent(out) :: croptype_frac(LDT_rc%lnc(n),LDT_rc%lnr(n))

! MIRCA2000 monthly_growing_area_grids crop area map (5 arcmin):
   integer, parameter :: IN_cols = 4320
   integer, parameter :: IN_rows = 2160
   real,    parameter :: IN_xres = 1.0/12.0 ! or 360deg/4320 col points
   real,    parameter :: IN_yres = 1.0/12.0 ! or 180deg/2160 row points
   real, parameter    :: PI = 3.1415926535898
   real, parameter    :: radius = 6371.0 ! [km]

   integer :: glpnc, glpnr             ! Parameter (global) total columns and rows
   integer :: subpnc, subpnr           ! Parameter subsetted columns and rows
   real    :: param_gridDesc(20)       ! Input parameter grid desc fgrd
   real    :: subparam_gridDesc(20)    ! Input parameter grid desc array
   integer, allocatable  :: lat_line(:,:), lon_line(:,:)

   logical :: file_exists
   integer :: ftn, ierr, nrec, ios 
   integer :: i, j, t, c, r, m, nc, nr
   integer :: mi                       ! Total number of input param grid array points
   integer :: mo                       ! Total number of output LIS grid array points
   real,    allocatable  :: gi1(:)     ! Input parameter 1d grid
   logical*1,allocatable :: li1(:)     ! Input logical mask (to match gi)
   real      :: go1(LDT_rc%lnc(n)*LDT_rc%lnr(n))            ! Output lis 1d grid
   logical*1 :: lo1(LDT_rc%lnc(n)*LDT_rc%lnr(n))            ! Output logical mask (to match go)

   real, allocatable :: var_read(:)           ! Read in parameter
   real, allocatable :: var_in(:,:,:)         ! Read in parameter
   real, allocatable :: read_crop_ann(:,:)    ! Read input annual mean

   real                                :: latc, lonc, area, d2r
!_______________________________________________________________________________________
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

!- Map Parameter Grid Info to LIS Target Grid/Projection/Subset Info:
   subparam_gridDesc = 0.
   call LDT_RunDomainPts( n, LDT_LSMCrop_struc(n)%crop_proj, &
        param_gridDesc, glpnc, glpnr, subpnc, subpnr, &
        subparam_gridDesc, lat_line, lon_line )


!- Perform check on grid and projection choices selected:
   if( LDT_rc%lis_map_proj == "latlon"  ) then 
     if( param_gridDesc(10) .gt. (LDT_rc%gridDesc(n,9)/LDT_rc%lis_map_resfactor)) then
        write(LDT_logunit,*) "[ERR] 'Native' MIRCA2000 crop files have been selected, "
        write(LDT_logunit,*) "    with a resolution (0.0833deg), but the LIS run domain resolution"
        write(LDT_logunit,*) "    selected is finer than that. Downscaling crop type information is "
        write(LDT_logunit,*) "    currently not supported in LDT."
        write(LDT_logunit,*) "Program stopping ..."
        call LDT_endrun
     endif
   else
    write(LDT_logunit,*) "[ERR] projection other than latlon is not supported "
    write(LDT_logunit,*) " due to only UpscaleByAveraging transformation "
    write(LDT_logunit,*) " is appropriate for MIRCA2000 crop map "
    write(LDT_logunit,*) "Program stopping ..."
    call LDT_endrun
   endif

   if( allocated(var_read) ) then
       deallocate(var_read)
   endif
   allocate( var_read(IN_cols*12) )
   allocate( var_in(subpnc,subpnr,12) )
!
!  READ IN MIRCA2000 CROP AREA FIELDS 
!
   ftn = LDT_getNextUnitNumber()
   open(ftn, file=trim(tempfile),status='old',form='unformatted',&
        access='direct',convert='little_endian',recl=4*IN_cols*12,iostat=ios)
!- Reverse-Y and read in 12-month values and subset
   nrec = 0
   do nr=subpnr, 1, -1
    nrec = IN_rows - lat_line(1,nr) + 1
    read(ftn,rec=nrec) var_read
    do nc = 1, subpnc
      do m = 1, 12
        var_in(nc,nr,m) = var_read((lon_line(nc,1)-1)*12 + m)
      end do
    end do
   end do

   deallocate( var_read )
   call LDT_releaseUnitNumber(ftn)

   write(LDT_logunit,*) "[INFO] Done reading MIRCA2000 crop type "

! -------------------------------------------------------------------
!     Compute annual maximum of the 12 month data
! -------------------------------------------------------------------
    allocate(read_crop_ann(subpnc,subpnr))
     do nr = 1, subpnr
       do nc = 1, subpnc
         read_crop_ann(nc,nr) = maxval(var_in(nc,nr,:))
         ! assign missing to non-existing crop type
         if( read_crop_ann(nc,nr) .le. 0. ) read_crop_ann(nc,nr) = LDT_rc%udef
       end do
     end do

! -------------------------------------------------------------------
!     AGGREGATING FINE-SCALE GRIDS TO COARSER LIS OUTPUT GRID
! -------------------------------------------------------------------
    mi = subpnc*subpnr
    allocate( gi1(mi), li1(mi) )
    gi1 = LDT_rc%udef;  li1 = .false.
    mo = LDT_rc%lnc(n)*LDT_rc%lnr(n)
    lo1 = .false.;  go1 = 0.

 !- Assign 2-D array to 1-D for aggregation routines:
    i = 0
    do nr = 1, subpnr
       do nc = 1, subpnc
          i = i + 1
          gi1(i) = read_crop_ann(nc,nr)
          if( gi1(i) < 0. ) gi1(i) = LDT_rc%udef
          if( gi1(i) .ne. LDT_rc%udef )  li1(i) = .true.
       enddo
    enddo
    deallocate( read_crop_ann )
    deallocate( var_in )

 !- Upscale routine:
    select case ( LDT_LSMCrop_struc(n)%crop_gridtransform )

   !- Single layer transformation:
      case( "none", "average", "neighbor", "bilinear", "budget-bilinear" )

       !- Transform parameter from original grid to LIS output grid:
          call LDT_transform_paramgrid(n, LDT_LSMCrop_struc(n)%crop_gridtransform, &
                   subparam_gridDesc, mi, 1, gi1, li1, mo, go1, lo1 )

       !- Convert 1D to 2D grid arrays and units from area (Ha) to fraction
          d2r = PI/180.
          j = 0
          do r = 1, LDT_rc%lnr(n)
             latc = (r-1) * LDT_rc%gridDesc(n,9) + LDT_rc%gridDesc(n,4)
             area = (sin(d2r*(latc+0.5*LDT_rc%gridDesc(n,9))) - sin(d2r*(latc-0.5*LDT_rc%gridDesc(n,9))))*(LDT_rc%gridDesc(n,10)*d2r)
             area = 100. * area * radius * radius ! in hectares as in MIRCA
             do c = 1, LDT_rc%lnc(n)
                j = j + 1
                if ( go1(j) < 0. ) then
                 croptype_frac(c,r) = 0.0
                else
                 croptype_frac(c,r) = go1(j)/area
                endif
              ! Make sure all negative values are set to universal undefined value:
                if( croptype_frac(c,r) < 0.  ) croptype_frac(c,r) = LDT_rc%udef
                if( croptype_frac(c,r) > 100.) croptype_frac(c,r) = LDT_rc%udef
             enddo
          enddo

     case default
       write(*,*)"[WARN] Other spatial grid transformations are not currently supported"
       write(*,*)"   for the tiled 'MIRCA2000' crop type maps. Please select either:"
       write(*,*)"  -- neighbor, bilinear, budget-bilinear (to downscale)"
       write(*,*)"  -- average (to upscale) or none"
       write(*,*)" Stopping program ..."
       call LDT_endrun

    end select  ! End vegtype/cnt aggregation method
    deallocate( gi1, li1 )
    deallocate( lat_line, lon_line )

 end subroutine readMIRCA2000Cropfiles

end module MIRCA2000_crops_module
