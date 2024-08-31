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
! !ROUTINE: read_VIC411_lc
!  \label{read_VIC411_lc}
!
! !REVISION HISTORY:
!  03 Sept 2004: Sujay Kumar; Initial Specification
!  23  Aug 2013: KR Arsenault; Implemented new features to read different
!                 resolutions and generate landmask from landcover
!
! !INTERFACE:
subroutine read_VIC411_lc(n, num_types, fgrd, maskarray)

! !USES:
  use LDT_coreMod,   only : LDT_rc
  use LDT_logMod,    only : LDT_logunit, LDT_getNextUnitNumber, &
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
!  This subroutine reads the VIC411 landcover data and returns the 
!  distribution of vegetation in each grid cell, in a lat/lon
!  projection.  Also, the landmask is either generated and/or 
!  read in this routine.
!
!  This VIC411-based land cover includes the following land cover classes:
!
!    layer0  -- water
!    layer1  -- evergreen needleleaf
!    layer2  -- evergreen broadleaf
!    layer3  -- deciduous needleleaf
!    layer4  -- deciduous broadleaf
!    layer5  -- mixed forest
!    layer6  -- woodland
!    layer7  -- wooded grassland
!    layer8  -- closed shrublands
!    layer9  -- open shrublands
!    layer10 -- grasslands
!    layer11 -- croplands
!    layer12 -- bare soil
!    (layer13) -- (should be urban fraction layer, but not used)
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
   integer, allocatable  :: lat_line(:,:)
   integer, allocatable  :: lon_line(:,:)
   real    :: isum, isum2

   real    :: vegperc(LDT_rc%lnc(n),LDT_rc%lnr(n),LDT_rc%nt)
!__________________________________________________________________

!- Assign additional land cover types, including generic water points: 
   select case ( LDT_rc%lc_type(n) )
    case ( "VIC", "UMD" )
      LDT_rc%bareclass    = 12
      LDT_rc%urbanclass   = 13
      LDT_rc%waterclass   = 14
      LDT_rc%snowclass    = 0
      LDT_rc%wetlandclass = 0
      LDT_rc%glacierclass = 0

    case default ! non-supported options
      write(LDT_logunit,*) "Land classification: ",trim(LDT_rc%lc_type(n)),&
                           " does not exist for VIC411 source."
      write(LDT_logunit,*) " -- Please select:  VIC411 "
      write(LDT_logunit,*) "program stopping ..."
      call LDT_endrun
   end select

   if( LDT_rc%mask_type(n) == "readin" ) then
     write(unit=LDT_logunit,fmt=*) " [INFO]  Currently for the VIC 4.1.1 UMD landcover "
     write(unit=LDT_logunit,fmt=*) "  map, 'create' is the only landmask option ... "
     write(unit=LDT_logunit,fmt=*) " Stopping ... "
     call LDT_endrun 
   endif
   vegperc  = 0.
   maskarray = 0.0
   fgrd     = 0.0

!- Check if land cover file exists:
   inquire( file=trim(LDT_rc%vfile(n)), exist=file_exists ) 
   if(.not. file_exists) then 
      write(LDT_logunit,*) "Landcover map: ",trim(LDT_rc%vfile(n))," does not exist "
      write(LDT_logunit,*) "program stopping ..."
      call LDT_endrun
   endif

!- Open LDT land cover file:
   write(unit=LDT_logunit,fmt=*)'[INFO] LAT/LON -- Reading ', trim(LDT_rc%vfile(n))
   ftn = LDT_getNextUnitNumber()
   open(ftn, file=LDT_rc%vfile(n), status='old', form='unformatted',&
        access ='direct', recl=4, iostat=ios1 )

! -------------------------------------------------------------------
!    PREPARE SUBSETTED PARAMETER GRID FOR READING IN NEEDED DATA
! -------------------------------------------------------------------
!- Map Parameter Grid Info to LIS Target Grid/Projection Info -- 
   subparam_gridDesc = 0.

   call LDT_RunDomainPts( n, LDT_rc%lc_proj, LDT_rc%lc_gridDesc(n,:), &
                  glpnc, glpnr, subpnc, subpnr,  &
                  subparam_gridDesc, lat_line, lon_line )

!- Determine if landcover map file is 2D or 3D:
   inquire( file=trim(LDT_rc%vfile(n)), size=file_size )
   if( file_size > glpnc*glpnr*4 ) then
      write(unit=LDT_logunit,fmt=*) '[INFO] Reading in a 3-D landcover map ...'
      file_dim = 3
   else
      file_dim = 2
   end if

!- Double-check tiled landcover if it is a tiled file:
   if( file_dim == 2 .and. LDT_rc%lc_gridtransform(n)=="tile" .and. &
       subparam_gridDesc(9) == (LDT_rc%gridDesc(n,9)/LDT_rc%lis_map_resfactor(n)) ) then
      write(LDT_logunit,*) " (in read_VIC411_lc) :: The 'tile' spatial transform option " 
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
         !- Dealing with UMD tiled classes -- current solution (KRA):
!            if( LDT_rc%lc_type(n) == "UMD" .and. t >= (LDT_rc%nt-1) ) exit  ! Global
            if( LDT_rc%lc_type(n) == "UMD" .and. t >= (LDT_rc%nt) ) exit    ! CONUS
            do r = 1, subpnr
               do c = 1, subpnc
                  line = (lat_line(c,r)-1)*glpnc + lon_line(c,r) + (t-1)*(glpnc*glpnr)
                  read(ftn,rec=line) vegperc(c,r,t)
                  if( vegperc(c,r,t) < 0. ) vegperc(c,r,t) = 0.
               enddo
            enddo
         enddo

      !- Locate land points and assign mask values:
         do r = 1, subpnr
            do c = 1, subpnc   
               isum = 0.
!               isum = sum( vegperc(c,r,1:LDT_rc%nt-2) )  ! Global
               isum = sum( vegperc(c,r,1:LDT_rc%nt-1) )   ! CONUS

               if( isum == 0. ) then         ! Ocean points
                 if( LDT_rc%mask_type(n) == "create" ) &
                    maskarray(c,r) = 0.
                  fgrd(c,r,LDT_rc%waterclass) = 1.0

               elseif( isum > 0.98 ) then   ! Landcover points

                  if( LDT_rc%mask_type(n) == "create" ) &
                    maskarray(c,r) = 1.
                    fgrd(c,r,LDT_rc%waterclass) = 0.
                    do t = 1, 12  
                       fgrd(c,r,t) = vegperc(c,r,t)
                    end do
                   
               end if

            enddo
         enddo

      !- "READ-IN" land mask file, if user-specified:
!         if( trim(LDT_rc%mask_type(n)) == "readin" ) then
!            if( LDT_rc%mask_source(n) == "VIC411" ) then
!               call read_VIC411_maskfile( n, fgrd, maskarray )
!            endif
!         endif

    ! ---
      else
         write(unit=LDT_logunit,fmt=*) " Resolutions of landcover file grid and LIS runtime grid "
         write(unit=LDT_logunit,fmt=*) "  are different and currently not supported ... "
         write(unit=LDT_logunit,fmt=*) " Stopping ... "
         call LDT_endrun 
      end if

    else
       write(LDT_logunit,*) ' No other land cover aggregation types are currently supported ... '
       write(LDT_logunit,*) ' -- program stopping ...'
       call LDT_endrun

    end if  ! End tiled option

    deallocate( lat_line, lon_line )

! -------------------------------------------------------------------
!    CREATE OR READ-IN (OR IMPOSE) LAND MASK FILE AND CREATE
!    SURFACE MAP
! -------------------------------------------------------------------

   call LDT_releaseUnitNumber(ftn)

end subroutine read_VIC411_lc
