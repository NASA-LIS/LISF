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
! !ROUTINE: read_ISA_lc
!  \label{read_ISA_lc}
!
! !REVISION HISTORY:
!  03 Sept 2004: Sujay Kumar; Initial Specification
!  23  May 2012: KR Arsenault; Implemented new features to read different
!                 resolutions and generate landmask from landcover
!
! !INTERFACE:
subroutine read_ISA_lc(n, num_types, fgrd, maskarray)

! !USES:
  use LDT_coreMod,     only : LDT_rc
  use LDT_logMod,      only : LDT_logunit, LDT_getNextUnitNumber, &
          LDT_releaseUnitNumber, LDT_verify, LDT_endrun
  use LDT_gridmappingMod    

  implicit none

! !ARGUMENTS: 
  integer, intent(in) :: n
  integer, intent(in) :: num_types
  real, intent(inout) :: fgrd(LDT_rc%lnc(n),LDT_rc%lnr(n),LDT_rc%nt)
  real, intent(inout) :: maskarray(LDT_rc%lnc(n),LDT_rc%lnr(n))
!
! !DESCRIPTION:
!  This subroutine reads the ISA landcover data and returns the 
!  distribution of vegetation in each grid cell, in a lat/lon
!  projection.  Also, the landmask is either generated and/or 
!  read in this routine.
!
!  This ISA-based land cover includes the following land cover classes:
!
!    layer0  -- water
!    layer1  -- evergreen broadleaf
!    layer2  -- deciduous broadleaf
!    layer3  -- mixed forest
!    layer4  -- evergreen needle
!    layer5  -- deciduous needle
!    layer6  -- savannas
!    layer7  -- grassland
!    layer8  -- urban and buildup
!    layer9  -- shrublands
!    layer10 -- tundra *******all zero because no tundra in USA******
!    layer11 -- barren
!    layer12 -- croplands
!
!  New:  Surface type array is generated here
!
!  The arguments are:
!  \begin{description}
!   \item[n]
!     index of nest
!   \item[fgrd]
!     fraction of grid covered by each vegetation type
!   \item[maskarray]
!     landmask for the region of interest
!   \end{description}
!EOP      

   integer :: ftn, ierr, ios1
   logical :: file_exists
   integer :: file_size
   integer :: file_dim
   integer :: i, t, c, r, line
   integer :: glpnc, glpnr             ! Parameter (global) total columns and rows
   integer :: subpnc, subpnr           ! Parameter subsetted columns and rows
   real    :: subparam_gridDesc(20)    ! Input parameter grid desc array
   integer, allocatable  :: lat_line(:,:), lon_line(:,:)
   integer :: isum, isum2
   real    :: readveg
   real    :: vegperc(LDT_rc%lnc(n),LDT_rc%lnr(n),LDT_rc%nt)
   real    :: vegtype(LDT_rc%lnc(n),LDT_rc%lnr(n))
!__________________________________________________________________

!- Assign additional land cover types, including generic water points: 
   select case ( LDT_rc%lc_type(n) )
    case ( "ISA" ) 
      LDT_rc%waterclass   = 13   ! Originally 0
      LDT_rc%urbanclass   = 8
      LDT_rc%bareclass    = 11
      LDT_rc%snowclass    = 0
      LDT_rc%wetlandclass = 0
      LDT_rc%glacierclass = 0
    case default ! non-supported options
      write(LDT_logunit,*) "Land classification: ",trim(LDT_rc%lc_type(n)),&
                           " does not exist for ISA source."
      write(LDT_logunit,*) " -- Please select:  ISA "
      write(LDT_logunit,*) "program stopping ..."
      call LDT_endrun
   end select

   vegperc  = 0.
   vegtype  = float(LDT_rc%waterclass)
   maskarray = 0.0
   fgrd     = 0.0

!- Check if land cover file exists:
   inquire( file=trim(LDT_rc%vfile(n)), exist=file_exists ) 
   if(.not. file_exists) then 
      write(LDT_logunit,*) "Landcover map: ",trim(LDT_rc%vfile(n))," does not exist "
      write(LDT_logunit,*) "program stopping ..."
      call LDT_endrun
   endif
   write(unit=LDT_logunit,fmt=*) "[INFO] Reading landcover file: ",trim(LDT_rc%vfile(n))

!- Open LDT land cover file:
   ftn = LDT_getNextUnitNumber()
   open(ftn, file=LDT_rc%vfile(n), status='old', form='unformatted',&
        access ='direct', recl=4, iostat=ios1, convert="little_endian")

! -------------------------------------------------------------------
!    PREPARE SUBSETTED PARAMETER GRID FOR READING IN NEEDED DATA
! -------------------------------------------------------------------
!- Map Parameter Grid Info to LIS Target Grid/Projection Info -- 
   subparam_gridDesc = 0.
   call LDT_RunDomainPts( n, LDT_rc%lc_proj, LDT_rc%lc_gridDesc(n,:), &
                  glpnc, glpnr, subpnc, subpnr,  &
                  subparam_gridDesc, lat_line, lon_line )

! -------------------------------------------------------------------

!- Only LAT/LON GCS currently supported:
   if( LDT_rc%lis_map_proj(n) .ne. "latlon" ) then
      write(LDT_logunit,*) "[INFO] For the ISA landcover type, the Lat/lon coordinate system"
      write(LDT_logunit,*) "     is only supported at this time.  Please select 'latlon' for both"
      write(LDT_logunit,*) "     LIS run domain and ISA landcover parameter domain." 
      write(LDT_logunit,*) " Stopping program ..."
      call LDT_endrun
   endif

!- Determine if landcover map file is 2D or 3D:
   inquire( file=trim(LDT_rc%vfile(n)), size=file_size )
   if( file_size > glpnc*glpnr*4 ) then
      write(LDT_logunit,*) "[INFO] Reading in a 3-D landcover map ..."
      file_dim = 3
   else
      file_dim = 2
   end if

!- Double-check tiled landcover if it is a tiled file:
   if( file_dim == 2 .and. LDT_rc%lc_gridtransform(n)=="tile" .and. &
       subparam_gridDesc(9) == (LDT_rc%gridDesc(n,9)/LDT_rc%lis_map_resfactor(n)) ) then
      write(LDT_logunit,*) " (in read_ISA_lc) :: The 'tile' spatial transform option " 
      write(LDT_logunit,*) "          has been selected, but the landcover file being read in"
      write(LDT_logunit,*) "          is not in vegetation tile-format order, and both your "
      write(LDT_logunit,*) "          landcover parameter and LIS run domain resolutions are "
      write(LDT_logunit,*) "          the same.  Please select option 'none' and run LDT again. "
      write(LDT_logunit,*) " Program stopping ..."
      call LDT_endrun
   end if

! -------------------------------------------------------------------
!    READ IN LAND COVER PARAMETER FIELDS (NON-TILED/TILED OPTIONS)
! -------------------------------------------------------------------

 !- (1) TILED land cover files:
   if ( LDT_rc%lc_gridtransform(n) == "tile" ) then

   !- Input parameter grid RES == LIS target grid RES:
      if( subparam_gridDesc(9) == (LDT_rc%gridDesc(n,9)/LDT_rc%lis_map_resfactor(n)) ) then

         line = 0
         do t = 1, LDT_rc%nt
            do r = 1, subpnr
               do c = 1, subpnc
                  line = (lat_line(c,r)-1)*glpnc + lon_line(c,r) + (t-1)*(glpnc*glpnr)
                  read(ftn,rec=line) readveg
    
               !- Reorder vegetation layers (1-12, 13=water)
                  if( t == 1 ) then     ! If water point, assign to 13
                     vegperc(c,r,13) = readveg
                  elseif( t > 1 ) then  ! For all other valid veg classes
                     vegperc(c,r,t-1) = readveg
                  end if

               enddo
            enddo
         enddo

      !- Locate land points and assign mask values:
         do r = 1, subpnr
            do c = 1, subpnc   
               isum = 0
               isum = sum( vegperc(c,r,1:LDT_rc%nt) ) 

               if( isum == 0 ) then         ! Ocean points
                 if( LDT_rc%mask_type(n) == "create" ) &
                    maskarray(c,r) = 0.
                  fgrd(c,r,LDT_rc%waterclass) = 1.0

               elseif( isum > 98 ) then   ! Landcover points

                ! Account for interior water points:
                  if( (vegperc(c,r,LDT_rc%waterclass)/100.) > LDT_rc%gridcell_water_frac(n)) then
                    if( LDT_rc%mask_type(n) == "create" ) &
                       maskarray(c,r) = 0.
                     fgrd(c,r,LDT_rc%waterclass) = 1.0
                     do t = 1, 12  
                        fgrd(c,r,t) = 0.
                     end do

                ! Account for land cover points:
                  else
                   ! 100% landcover percentage:
                    if( LDT_rc%mask_type(n) == "create" ) &
                      maskarray(c,r) = 1.
                     fgrd(c,r,LDT_rc%waterclass) = 0.

                     isum2 = 0
                     isum2 = sum( vegperc(c,r,1:LDT_rc%nt-1) ) 
                     do t = 1, 12  
                        fgrd(c,r,t) = vegperc(c,r,t)/float(isum2)
                     end do
                  end if

               end if
            enddo
         enddo

      !- "READ-IN" land mask file, if user-specified:
         if( trim(LDT_rc%mask_type(n)) == "readin" ) then
            if( LDT_rc%mask_source(n) == "ISA" ) then
               call read_ISA_maskfile( n, fgrd, maskarray )
            else
               call read_maskfile( n, vegtype, fgrd, maskarray )
            endif
         endif

    ! ---
      else
         write(unit=LDT_logunit,fmt=*) "[INFO] Resolutions of ISA landcover file grid and LIS runtime grid "
         write(unit=LDT_logunit,fmt=*) "  are different and different resolutions are currently not supported ... "
         write(unit=LDT_logunit,fmt=*) " Program stopping ... "
         call LDT_endrun 
      end if

    else
       write(LDT_logunit,*) "[INFO] No other landcover spatial transformation methods are available "
       write(LDT_logunit,*) "  at this time for the ISA land classification type ... "
       write(LDT_logunit,*) " Program stopping ..."
       call LDT_endrun

    end if  ! End tiled option
    deallocate( lat_line, lon_line )

! -------------------------------------------------------------------
!    CREATE OR READ-IN (OR IMPOSE) LAND MASK FILE AND CREATE
!    SURFACE MAP
! -------------------------------------------------------------------

   call LDT_releaseUnitNumber(ftn)

end subroutine read_ISA_lc
