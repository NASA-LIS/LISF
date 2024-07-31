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
! !ROUTINE: read_UMDCROPMAP_croptype
!  \label{read_UMDCROPMAP_croptype}
!
! !REVISION HISTORY:
!  03 Sept 2004: Sujay Kumar; Initial Specification
!  23  May 2012: KR Arsenault; Implemented new features to read different
!                 resolutions and generate landmask from landcover
!
! !INTERFACE:
subroutine read_UMDCROPMAP_croptype(n, num_types, fgrd)

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
!  This subroutine reads the UMD+CROPMAP landcover data and returns the 
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

   integer :: ftn, ierr, ios1
   logical :: file_exists
   integer :: i, t, c, r, line
   integer :: glpnc, glpnr             ! Parameter (global) total columns and rows
   integer :: subpnc, subpnr           ! Parameter subsetted columns and rows
   integer :: mi                       ! Total number of input param grid array points
   integer :: mo                       ! Total number of output LIS grid array points
   real    :: param_gridDesc(20)       ! Input parameter grid desc fgrd
   real    :: subparam_gridDesc(20)    ! Input parameter grid desc array
   integer, allocatable  :: lat_line(:,:), lon_line(:,:)
   real,    allocatable  :: gi(:)      ! Input parameter 1d grid
   logical*1,allocatable :: li(:)      ! Input logical mask (to match gi)

   real      :: go1(LDT_rc%lnc(n)*LDT_rc%lnr(n))            ! Output lis 1d grid
   logical*1 :: lo1(LDT_rc%lnc(n)*LDT_rc%lnr(n))            ! Output logical mask (to match go)
   real      :: go2(LDT_rc%lnc(n)*LDT_rc%lnr(n), num_types) ! Output lis 1d grid
   logical*1 :: lo2(LDT_rc%lnc(n)*LDT_rc%lnr(n), num_types) ! Output logical mask (to match go)

   integer   :: file_size
   integer   :: file_dim
   integer   :: water_class
   real      :: crop_array(num_types)
   real      :: vegcnt(LDT_rc%lnc(n),LDT_rc%lnr(n),num_types)
   real      :: vegtype(LDT_rc%lnc(n),LDT_rc%lnr(n))

!__________________________________________________________________

   water_class  = 33

   vegcnt  = 0.
   vegtype = float(water_class)
   crop_array = 0.0
   fgrd    = 0.0

!- Set parameter grid array inputs:
   LDT_LSMCrop_struc(n)%crop_proj = "latlon"
   param_gridDesc(1)  = 0.          ! Latlon
   param_gridDesc(2)  = 464
   param_gridDesc(3)  = 224
   param_gridDesc(4)  = 25.0625     ! LL lat 
   param_gridDesc(5)  = -124.9375   ! LL lon 
   param_gridDesc(6)  = 128
   param_gridDesc(7)  = 52.9375     ! UR lat
   param_gridDesc(8)  = -67.0625    ! UR lon
   param_gridDesc(9)  = 0.125
   param_gridDesc(10) = 0.125
   param_gridDesc(20) = 64

!- Check if land cover file exists:
   inquire( file=trim(LDT_LSMCrop_struc(n)%croptfile), exist=file_exists ) 
   if(.not. file_exists) then 
      write(LDT_logunit,*) "Crop type map: ",trim(LDT_LSMCrop_struc(n)%croptfile)," does not exist. "
      write(LDT_logunit,*) "program stopping ..."
      call LDT_endrun
   endif
   write(LDT_logunit,*)"[INFO] Reading UMD/CROPMAP: ",trim(LDT_LSMCrop_struc(n)%croptfile)
   write(LDT_logunit,*)"** NOTE: Only 'mode' spatial transform currently works"
   write(LDT_logunit,*)"         with the 'UMDCROPMAP' crop type map source." 

!- Open LDT land cover file:
   ftn = LDT_getNextUnitNumber()
   open(ftn, file=LDT_LSMCrop_struc(n)%croptfile, status='old', form='unformatted',&
        access ='direct', recl=4, iostat=ios1)

! -------------------------------------------------------------------
!    PREPARE SUBSETTED PARAMETER GRID FOR READING IN NEEDED DATA
! -------------------------------------------------------------------
!- Map Parameter Grid Info to LIS Target Grid/Projection Info -- 
   subparam_gridDesc = 0.
   call LDT_RunDomainPts( n, LDT_LSMCrop_struc(n)%crop_proj, param_gridDesc, &
                  glpnc, glpnr, subpnc, subpnr,  &
                  subparam_gridDesc, lat_line, lon_line )

!- Determine if crop map file is 2D or 3D:
   inquire( file=trim(LDT_LSMCrop_struc(n)%croptfile), size=file_size )
   if( file_size > glpnc*glpnr*4 ) then
      write(unit=LDT_logunit,fmt=*) '[INFO] Opening/Reading 3-D crop map ...'
      file_dim = 3   ! 3-D fields
   else 
      file_dim = 2   ! 2-D fields
   end if

   if( file_dim == 3 .and. LDT_LSMCrop_struc(n)%crop_gridtransform == "mode" ) then  
     if( subparam_gridDesc(9) .ne. (LDT_rc%gridDesc(n,9)/LDT_rc%lis_map_resfactor(n)) ) then
       write(*,*) "[WARN] IF 3-D TILED CROP MAP DOES NOT HAVE THE SAME RESOLUTION "
       write(*,*) "  AS THE LIS RUN-DOMAIN, THEN 'MODE' OPTION CANNOT BE SELECTED,"
       write(*,*) "  OR YOU NEED TO SELECT 'TILE' AS YOUR OPTION."
       write(*,*) " Stopping ..."
       call LDT_endrun
     endif
   endif

! -------------------------------------------------------------------
!    READ IN LAND COVER PARAMETER FIELDS (NON-TILED/TILED OPTIONS)
! -------------------------------------------------------------------

!- Initialize parameter read-in array:
  select case ( LDT_LSMCrop_struc(n)%crop_gridtransform ) 

!    case ( "mode", "tile" )
    case ( "mode" )
      line = 0
      do t = 1, num_types
         do r = 1, subpnr
            do c = 1, subpnc
               line = (lat_line(c,r)-1)*glpnc + lon_line(c,r) + (t-1)*(glpnc*glpnr)
               read(ftn,rec=line) vegcnt(c,r,t)
            enddo
         enddo
      enddo

    case default
       write(LDT_logunit,*)" No other aggregation types are currently supported for UMD/CROPLAND." 
       write(LDT_logunit,*)" -- Program stopping ..."
       call LDT_endrun
    end select  
    deallocate( lat_line, lon_line )

! -------------------------------------------------------------------
!     AGGREGATING FINE-SCALE GRIDS TO COARSER LIS OUTPUT GRID
! -------------------------------------------------------------------
   mi = subpnc*subpnr
   mo = LDT_rc%lnc(n)*LDT_rc%lnr(n)
   allocate( gi(mi), li(mi) )
   gi = float(water_class)
   li = .false.
   lo1 = .false.;  lo2 = .false.

   select case( LDT_LSMCrop_struc(n)%crop_gridtransform ) 

  !- (a) Estimate dominant landcover/crop types:
     case( "mode" )

    !- Aggregate:
       if( file_dim == 2 ) then   ! 2-D fields
  
       !- Transform parameter from original grid to LIS output grid:
          call LDT_transform_paramgrid(n, LDT_LSMCrop_struc(n)%crop_gridtransform, &
                   subparam_gridDesc, mi, 1, gi, li, mo, go1, lo1 )

       !- Convert 1D dominant veg type to 2D grid arrays:
          i = 0
          do r = 1, LDT_rc%lnr(n)
             do c = 1, LDT_rc%lnc(n)
                i = i + 1
                vegtype(c,r) = go1(i)
             enddo
          enddo
 
    !- No spatial aggregation needed; just dom. class selected:
       elseif( file_dim == 3 ) then   ! 3-D fields

          do r = 1, subpnr
             do c = 1, subpnc
                do t = 1, num_types
                   if( t > 13 ) crop_array(t) = vegcnt(c,r,t)
                end do
                if( sum(crop_array(14:num_types)) == 0. ) then
                  fgrd(c,r,1) = LDT_rc%udef
                else
                  fgrd(c,r,1) = maxloc(crop_array,1)
                end if
             enddo
          enddo
          return  ! Return dominant vegtype locations as a class to core routine

       endif  ! End 2D/3D check

!     case( "tile" )
     case default
       write(*,*) "[WARN] Other spatial grid transformations are not currently supported  "
       write(*,*) "  for the tiled UMD-CROPLAND landcover/crop map type.  Please select either:"
!       write(*,*) "  -- Mode or Tile "
       write(*,*) "  -- Mode "
       write(*,*) " Stopping ..."
       call LDT_endrun
   end select  ! End vegtype/cnt aggregation method
   deallocate( gi, li )

! ........................................................................

!- Bring 2-D Vegtype to 3-D Vegcnt tile space:
   if ( LDT_LSMCrop_struc(n)%crop_gridtransform == "none" .or. &
        LDT_LSMCrop_struc(n)%crop_gridtransform == "mode" ) then  ! -- NON-TILED SURFACES
     do r = 1, LDT_rc%lnr(n)
        do c = 1, LDT_rc%lnc(n)
           if ( vegtype(c,r) .le. 0 ) then
              vegtype(c,r) = float(water_class)
           endif
           if ( (nint(vegtype(c,r)) .ne. water_class ) .and. &
                (nint(vegtype(c,r)) .ne. LDT_rc%udef)) then
              vegcnt(c,r,NINT(vegtype(c,r))) = 1.0
           endif
        enddo
     end do
   endif   ! End NON-TILED vegetation option

!- Estimate fraction of grid (fgrid) represented by vegetation type::
   call param_index_fgrdcalc( n, LDT_LSMCrop_struc(n)%crop_proj, &
        LDT_LSMCrop_struc(n)%crop_gridtransform, &
        water_class, num_types, vegcnt, fgrd )

   call LDT_releaseUnitNumber(ftn)

   write(LDT_logunit,*) "[INFO] Done reading UMD/CROPMAP. "


end subroutine read_UMDCROPMAP_croptype
