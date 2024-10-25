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
! !ROUTINE: read_USGS_lc
!  \label{read_USGS_lc}
!
! !REVISION HISTORY:
!  03 Sept 2004: Sujay Kumar; Initial Specification
!  23  May 2012: KR Arsenault; Implemented new features to read different
!                 resolutions and generate landmask from landcover
!
! !INTERFACE:
subroutine read_USGS_lc(n, num_types, fgrd, maskarray)

! !USES:
  use LDT_coreMod,  only : LDT_rc
  use LDT_logMod,   only : LDT_logunit, LDT_getNextUnitNumber, &
       LDT_releaseUnitNumber, LDT_verify, LDT_endrun
  use LDT_gridmappingMod    
  use LDT_fileIOMod
  use LDT_paramTileInputMod, only: param_index_fgrdcalc

  implicit none

! !ARGUMENTS: 
  integer, intent(in) :: n
  integer, intent(in) :: num_types
  real, intent(inout) :: fgrd(LDT_rc%lnc(n),LDT_rc%lnr(n),LDT_rc%nt)
  real, intent(inout) :: maskarray(LDT_rc%lnc(n),LDT_rc%lnr(n))
!
! !DESCRIPTION:
!  This subroutine reads the USGS landcover data and returns the 
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
!   \item[maskarray]
!     landmask for the region of interest
!   \end{description}
!EOP      
!
   integer :: ftn, ierr, ios1
   logical :: file_exists
   integer :: file_size
   integer :: file_dim
   integer :: i, t, c, r, k, line
   integer :: nt_adj
   integer :: glpnc, glpnr             ! Parameter (global) total columns and rows
   integer :: subpnc, subpnr           ! Parameter subsetted columns and rows
   integer :: mi                       ! Total number of input param grid array points
   integer :: mo                       ! Total number of output LIS grid array points
   real    :: subparam_gridDesc(20)    ! Input parameter grid desc array
   integer, allocatable  :: lat_line(:,:), lon_line(:,:)
   real,    allocatable  :: read_inputparm(:,:)  ! Read input parameter
   real,    allocatable  :: read_vegcnt(:,:,:)   ! Read input parameter
   real,    allocatable  :: gi(:)      ! Input parameter 1d grid
   logical*1,allocatable :: li(:)      ! Input logical mask (to match gi)

   real      :: go1(LDT_rc%lnc(n)*LDT_rc%lnr(n))            ! Output lis 1d grid
   logical*1 :: lo1(LDT_rc%lnc(n)*LDT_rc%lnr(n))            ! Output logical mask (to match go)
   real      :: go2(LDT_rc%lnc(n)*LDT_rc%lnr(n), LDT_rc%nt) ! Output lis 1d grid
   logical*1 :: lo2(LDT_rc%lnc(n)*LDT_rc%lnr(n), LDT_rc%nt) ! Output logical mask (to match go)
   real      :: vegcnt(LDT_rc%lnc(n),LDT_rc%lnr(n),LDT_rc%nt)
   real      :: vegtype(LDT_rc%lnc(n),LDT_rc%lnr(n))
!__________________________________________________________________

!- Assign additional land cover types, including generic water points: 
   select case ( LDT_rc%lc_type(n) )
    case ( "USGS" )
      LDT_rc%urbanclass   = 1
      LDT_rc%waterclass   = 16
      LDT_rc%bareclass    = 19
      LDT_rc%snowclass    = 24
      LDT_rc%wetlandclass = 25
      LDT_rc%glacierclass = 26

    case default ! non-supported options
      write(LDT_logunit,*) "The land classification: ",trim(LDT_rc%lc_type(n)),&
                           " does not exist for USGS source."
      write(LDT_logunit,*) " -- Please select:  USGS "
      write(LDT_logunit,*) "Program stopping ..."
      call LDT_endrun
   end select

!- Check if land cover file exists:
   inquire( file=trim(LDT_rc%vfile(n)), exist=file_exists ) 
   if(.not. file_exists) then 
      write(LDT_logunit,*) "Landcover map: ",trim(LDT_rc%vfile(n))," does not exist. "
      write(LDT_logunit,*) "Program stopping ..."
      call LDT_endrun
   endif
   write(unit=LDT_logunit,fmt=*)"[INFO] Reading USGS (LIS team-derived) &
                                 landcover file: ", trim(LDT_rc%vfile(n))

!- Open LDT land cover file:
   ftn = LDT_getNextUnitNumber()
   open(ftn, file=LDT_rc%vfile(n), status='old', form='unformatted',&
       access ='direct', recl=4, iostat=ios1)

! -------------------------------------------------------------------

   vegcnt   = 0.
   vegtype  = float(LDT_rc%waterclass)
   maskarray = 0.0
   fgrd     = 0.0

! -------------------------------------------------------------------
!- Map Parameter Grid Info to LIS Target Grid/Projection Info -- 
   subparam_gridDesc = 0.
   call LDT_RunDomainPts( n, LDT_rc%lc_proj, LDT_rc%lc_gridDesc(n,:), &
                          glpnc, glpnr, subpnc, subpnr,  &
                          subparam_gridDesc, lat_line, lon_line )

! -------------------------------------------------------------------

!- Determine if landcover map file is 2D or 3D:
   nt_adj = LDT_rc%nt

   inquire( file=trim(LDT_rc%vfile(n)), size=file_size )
   if( file_size >= glpnc*glpnr*nt_adj*4 ) then
      write(unit=LDT_logunit,fmt=*) "[INFO] Reading a 3-D landcover map ..."
      file_dim = 3
   else
      write(unit=LDT_logunit,fmt=*) "[INFO] Reading a 2-D landcover map ..."
      file_dim = 2
   end if

!- Exceptions for Latlon/Lambert/Mercator projections:
   if( LDT_rc%lis_map_proj(n) == "latlon"   .or. &
       LDT_rc%lis_map_proj(n) == "mercator" .or. &
       LDT_rc%lis_map_proj(n) == "lambert" ) then

  !- File-size check for 2KM and 3KM global files:
     if( (subparam_gridDesc(9) < 0.04  .and. &
          subparam_gridDesc(9) > 0.01 ) .or. &
         (subparam_gridDesc(10) < 0.04 .and. &
          subparam_gridDesc(10) > 0.01) ) then
        if ( file_size < 0. )  file_dim = 3
     endif
  !- Double-check tiled landcover if it is a tiled file:
     if( file_dim == 2 .and. LDT_rc%lc_gridtransform(n)=="tile" ) then
       if( subparam_gridDesc(9) == (LDT_rc%gridDesc(n,9)/LDT_rc%lis_map_resfactor(n)) .or. &
           subparam_gridDesc(10) == (LDT_rc%gridDesc(n,10)/LDT_rc%lis_map_resfactor(n))) then
        write(LDT_logunit,*) "[ERR] in read_USGS_lc :: The 'tile' spatial transform option " 
        write(LDT_logunit,*) "   has been selected, but the landcover file being read in"
        write(LDT_logunit,*) "   is not in vegetation tile-format order, and both your "
        write(LDT_logunit,*) "   landcover parameter and LIS run domain resolutions are "
        write(LDT_logunit,*) "   the same.  Please select option 'none' and run LDT again. "
        write(LDT_logunit,*) " Program stopping ..."
        call LDT_endrun
      end if
     endif
  !- Lat-lon Landcover files cannot have "none" for spatial transform if
  !   LIS run domain projection is different:
     if((LDT_rc%lis_map_proj(n) == "mercator" .or.  &
         LDT_rc%lis_map_proj(n) == "lambert") .and. &
         LDT_rc%lc_gridtransform(n) == "none" ) then
        write(LDT_logunit,*) "[ERR] Landcover file has lat-lon grid coordinate system and"
        write(LDT_logunit,*) "  being translated to different output grid and projection."
        write(LDT_logunit,*) "  Therefore, the spatial transform option 'none' cannot be used here."
        write(LDT_logunit,*) "  Please then select another option (e.g., 'neighbor' for when input and "
        write(LDT_logunit,*) "  output grid resolutions are similar. "
        write(LDT_logunit,*) "Program stopping ..."
        call LDT_endrun
     endif

   endif   ! End other LIS-grid projection checks

! -------------------------------------------------------------------
!    READ IN LAND COVER PARAMETER FIELDS (NON-TILED/TILED OPTIONS)
! -------------------------------------------------------------------

!- Initialize parameter read-in array:
   allocate( read_inputparm(subpnc, subpnr) )
   read_inputparm = float(LDT_rc%waterclass)
   line = 0

   select case( LDT_rc%lc_gridtransform(n) ) 

     case( "none", "neighbor", "mode", "tile" ) 

 !- (1) TILED land cover files where INPUT Resolution MATCHES LIS Domain Resolution:
 !     (Old LIS-6 handling)
      if( file_dim == 3 ) then   ! 3-D file-size

         allocate( read_vegcnt(subpnc,subpnr,LDT_rc%nt) )
         read_vegcnt = 0.

         do t = 1, LDT_rc%nt
            do r = 1, subpnr
               do c = 1, subpnc
                  line = (lat_line(c,r)-1)*glpnc + lon_line(c,r) + (t-1)*(glpnc*glpnr)
                  read(ftn,rec=line) read_vegcnt(c,r,t)
               enddo
            enddo
         enddo

      !- Map readin 3-D vegcnt to output 3-D vegcnt:
         if( LDT_rc%lis_map_proj(n) == "latlon" ) then  ! Latlon only
           vegcnt = read_vegcnt

       ! Other projections ....
         else    ! Lambert, Gaussian, etc.
           write(LDT_logunit,*) "[INFO] It is NOT recommended to use the LIS-based 3D tiled vegetation"
           write(LDT_logunit,*) "     files for map projections (e.g., Lambert, Gaussian, etc.)"
           write(LDT_logunit,*) "     other than Lat/Lon."

         ! Use nearest neighbor for Lambert or other projections to read in 3D tiled veg.
           call veg3d_neighbor( n, subparam_gridDesc, subpnc, subpnr, read_vegcnt, vegcnt  )

         endif  ! End LDT run domain projection check

         deallocate( read_vegcnt )

   
    !- (2) Reading in 2D landcover map:
       else   ! 2-D
         do r = 1, subpnr
            do c = 1, subpnc
               line = (lat_line(c,r)-1)*glpnc + lon_line(c,r)
               read(ftn,rec=line) read_inputparm(c,r)

            !- Reassign 0 water class to an assigned classification value:
               if( read_inputparm(c,r) == 0. )  &
                   read_inputparm(c,r) = float(LDT_rc%waterclass)

            enddo
         enddo
       endif 

  !- Case when incorrect/unrecognized spatial transformation is selected:
     case default
       write(LDT_logunit,*)" No other land cover aggregation types are currently supported ... "
       write(LDT_logunit,*)" -- Program stopping ..."
       call LDT_endrun
    end select 
    deallocate( lat_line, lon_line )

! -------------------------------------------------------------------
!     AGGREGATING FINE-SCALE GRIDS TO COARSER LIS OUTPUT GRID
! -------------------------------------------------------------------

 !- No spatial transform applied:
    if( LDT_rc%lc_gridtransform(n) == "none" .and. &
        LDT_rc%lis_map_proj(n) == "latlon" ) then

       write(LDT_logunit,*) " No aggregation applied for parameter file ... "
       vegtype(:,:) = read_inputparm(:,:)

 !- Other transformations applied:
    else

       mi = subpnc*subpnr
       mo = LDT_rc%lnc(n)*LDT_rc%lnr(n)
       allocate( gi(mi), li(mi) )
       gi = float(LDT_rc%waterclass)
       li = .false.; lo1 = .false.;  lo2 = .false.

    !- Assign 2-D array to 1-D for aggregation routines:
       i = 0
       do r = 1, subpnr
          do c = 1, subpnc;  i = i + 1
             gi(i) = read_inputparm(c,r)
             if( gi(i) .ne. LDT_rc%udef ) li(i) = .true.
          enddo
       enddo

    !- (a) Estimate NON-TILED dominant land cover types (vegtype):
       if ( LDT_rc%lc_gridtransform(n) == "mode" .or. &
            LDT_rc%lc_gridtransform(n) == "neighbor" ) then

       !- Transform parameter from original grid to LIS output grid:
          call LDT_transform_paramgrid(n, LDT_rc%lc_gridtransform(n), &
                   subparam_gridDesc, mi, 1, gi, li, mo, go1, lo1 )

       !- Convert 1D dominant veg type to 2D grid arrays:
          i = 0
          do r = 1, LDT_rc%lnr(n)
             do c = 1, LDT_rc%lnc(n)
                i = i + 1
                vegtype(c,r) = go1(i)
             enddo
          enddo

    !- (b) Estimate TILED land cover files (vegcnt):
       elseif( LDT_rc%lc_gridtransform(n) == "tile" .and. &
               file_dim == 2 ) then   ! 2-d map

       !- Transform parameter from original grid to LIS output grid:
          call LDT_transform_paramgrid(n, LDT_rc%lc_gridtransform(n), &
                   subparam_gridDesc, mi, LDT_rc%nt, gi, li, mo, go2, lo2 )

       !- Convert 1D vegcnt to 2D grid arrays:
          i = 0
          do r = 1, LDT_rc%lnr(n) 
             do c = 1, LDT_rc%lnc(n)  
                i = i + 1
                do t = 1, LDT_rc%nt
                   vegcnt(c,r,t) = go2(i,t)
                end do
             enddo
          enddo

       endif
       deallocate( gi, li )

    endif  ! End vegtype/cnt aggregation method

    deallocate( read_inputparm )

! ........................................................................

!- Bring 2-D Vegtype to 3-D Vegcnt tile space:
   if ( LDT_rc%lc_gridtransform(n) == "none" .or. &
        LDT_rc%lc_gridtransform(n) == "mode" .or. &  
        LDT_rc%lc_gridtransform(n) == "neighbor" ) then  ! -- NON-TILED SURFACES
     do r = 1, LDT_rc%lnr(n)
        do c = 1, LDT_rc%lnc(n)
           if ( vegtype(c,r) .le. 0 ) then
              vegtype(c,r) = float(LDT_rc%waterclass)
           endif
           if ( (nint(vegtype(c,r)) .ne. LDT_rc%waterclass ) .and. &
                (nint(vegtype(c,r)) .ne. LDT_rc%udef)) then
              vegcnt(c,r,NINT(vegtype(c,r))) = 1.0
           endif
        enddo
     end do
   endif   ! End NON-TILED vegetation option

!- Estimate fraction of grid (fgrid) represented by vegetation type::

   call param_index_fgrdcalc( n, LDT_rc%lc_proj, LDT_rc%lc_gridtransform(n), &
                              LDT_rc%waterclass, LDT_rc%nt, vegcnt, fgrd )


! -------------------------------------------------------------------
!    CREATE OR READ-IN (OR IMPOSE) LAND MASK FILE AND CREATE
!    SURFACE MAP
! -------------------------------------------------------------------

!- "READ-IN" land mask file, if user-specified:
   if( LDT_rc%mask_type(n) == "readin" ) then

     call read_maskfile( n, vegtype, fgrd, maskarray )


!- "CREATE" land mask and surface type fields (user-specified):
   elseif( LDT_rc%mask_type(n) == "create" ) then

      call create_maskfile( n, LDT_rc%nt, LDT_rc%lc_gridtransform(n), &
                  vegtype, vegcnt, maskarray )

   end if

   call LDT_releaseUnitNumber(ftn)

end subroutine read_USGS_lc

