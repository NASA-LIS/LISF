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
! !ROUTINE: read_USGSNative_lc
!  \label{read_USGSNative_lc}
!
! !REVISION HISTORY:
!  03 Sept 2004: Sujay Kumar; Initial Specification
!  23  May 2012: KR Arsenault; Implemented new features to read different
!                 resolutions and generate landmask from landcover
!  30  Apr 2014: KR Arsenault; Added reader for "Native" USGS landcover
!
! !INTERFACE:
subroutine read_USGSNative_lc(n, num_types, fgrd, maskarray)

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
  integer, parameter :: input_cols = 360*120     ! 43200
  integer, parameter :: input_rows = 180*120     ! 21600
  real,    parameter :: input_xres = 1.0/120.0
  real,    parameter :: input_yres = 1.0/120.0
 ! 120 represents the number of gridcells per 1-deg gridbox

  integer :: ftn, ierr, ios1, nrec
  logical :: file_exists
  integer :: i, t, c, r, k, line
  integer :: glpnc, glpnr             ! Parameter (global) total columns and rows
  integer :: subpnc, subpnr           ! Parameter subsetted columns and rows
  integer :: mi                       ! Total number of input param grid array points
  integer :: mo                       ! Total number of output LIS grid array points
  real    :: param_gridDesc(20)       ! Input parameter grid desc fgrd
  real    :: subparam_gridDesc(20)    ! Input parameter for subsetted grid desc array
  integer, allocatable  :: lat_line(:,:), lon_line(:,:)
  real,    allocatable  :: gi(:)      ! Input parameter 1d grid
  logical*1,allocatable :: li(:)      ! Input logical mask (to match gi)

  real      :: go1(LDT_rc%lnc(n)*LDT_rc%lnr(n))            ! Output lis 1d grid
  logical*1 :: lo1(LDT_rc%lnc(n)*LDT_rc%lnr(n))            ! Output logical mask (to match go)
  real      :: go2(LDT_rc%lnc(n)*LDT_rc%lnr(n), LDT_rc%nt) ! Output lis 1d grid
  logical*1 :: lo2(LDT_rc%lnc(n)*LDT_rc%lnr(n), LDT_rc%nt) ! Output logical mask (to match go)
  real      :: vegcnt(LDT_rc%lnc(n),LDT_rc%lnr(n),LDT_rc%nt)
  real      :: vegtype(LDT_rc%lnc(n),LDT_rc%lnr(n))

  integer*1,allocatable  :: veg_int(:)
  real,     allocatable  :: read_inputparm(:,:)  ! Read input parameter
! __________________________________________________________________

!- Assign additional land cover types, including generic water points: 
   select case ( LDT_rc%lc_type(n) )
    case ( "USGS" )
      LDT_rc%urbanclass   = 1
      LDT_rc%waterclass   = 16
      LDT_rc%bareclass    = 19
      LDT_rc%snowclass    = 24
      LDT_rc%wetlandclass = 25
      LDT_rc%glacierclass = 26

   !- Set parameter grid array inputs:
      param_gridDesc(1)  = 0.             ! Latlon
      param_gridDesc(2)  = input_cols
      param_gridDesc(3)  = input_rows
      param_gridDesc(4)  = -90.0  + (input_yres/2) ! LL lat
      param_gridDesc(5)  = -180.0 + (input_xres/2) ! LL lon
      param_gridDesc(6)  = 128
      param_gridDesc(7)  =  90.0 - (input_yres/2)  ! UR lat
      param_gridDesc(8)  = 180.0 - (input_xres/2)  ! UR lon
      param_gridDesc(9)  = input_yres     ! dy: 0.0083333
      param_gridDesc(10) = input_xres     ! dx: 0.0083333
      param_gridDesc(20) = 64

    case default ! non-supported options
      write(LDT_logunit,*) "[ERR] The land classification: ",trim(LDT_rc%lc_type(n)),&
                           " does not exist for USGS source."
      write(LDT_logunit,*) " -- Please select:  USGS "
      write(LDT_logunit,*) "Program stopping ..."
      call LDT_endrun
   end select

!- Check if land cover file exists:
   inquire( file=trim(LDT_rc%vfile(n)), exist=file_exists ) 
   if(.not. file_exists) then 
      write(LDT_logunit,*) "[ERR] Landcover map: ",trim(LDT_rc%vfile(n))," does not exist. "
      write(LDT_logunit,*) "Program stopping ..."
      call LDT_endrun
   endif
   write(unit=LDT_logunit,fmt=*)"[INFO] Reading USGS-Native landcover file: ", trim(LDT_rc%vfile(n))

!- Open 30-arcsec USGS vegetation/land cover file:
   ftn = LDT_getNextUnitNumber()
   open(ftn, file=LDT_rc%vfile(n), status='old', form='unformatted',&
        access ='direct', recl=input_cols, iostat=ios1)  ! input_cols longitude

! -------------------------------------------------------------------

   vegtype   = float(LDT_rc%waterclass)
   vegcnt    = 0.
   fgrd      = 0.0
   maskarray = 0.0

! -------------------------------------------------------------------
!- Map Parameter Grid Info to LIS Target Grid/Projection Info -- 
   subparam_gridDesc = 0.
   call LDT_RunDomainPts( n, LDT_rc%lc_proj, param_gridDesc(:), &
            glpnc, glpnr, subpnc, subpnr, subparam_gridDesc, lat_line, lon_line )

! -------------------------------------------------------------------
!    READ IN LAND COVER PARAMETER FIELDS (NON-TILED/TILED OPTIONS)
! -------------------------------------------------------------------

!- Initialize parameter read-in array:
   allocate( read_inputparm(subpnc, subpnr) )
   read_inputparm = float(LDT_rc%waterclass)

   select case( LDT_rc%lc_gridtransform(n) ) 

     case( "none", "neighbor", "mode", "tile" ) 

    !- Reverse-Y and read in values:
       allocate( veg_int(input_cols) )
       nrec = 0
       do r = subpnr, 1, -1
          nrec = input_rows - lat_line(1,r) + 1
          read(ftn,rec=nrec,iostat=ierr) veg_int
          if(ierr /= 0) then
            write(*,*) "ERR: USGS-Native Landcover reader:: "
            write(*,*) "    Read error on record nrec = ", nrec
            write(*,*) " Stopping - read error ..."
            call LDT_endrun
          endif
          do c = 1, subpnc
             read_inputparm(c,r) = veg_int(lon_line(c,1))

          !- Reassign 0 water class to an assigned classification value:
             if( read_inputparm(c,r) == 0. ) then
                 read_inputparm(c,r) = float(LDT_rc%waterclass)
             endif
           enddo
        enddo
        deallocate( veg_int )

  !- Case when incorrect/unrecognized spatial transformation is selected:
     case default
       write(LDT_logunit,*) " No other land cover aggregation types are currently supported ... "
       write(LDT_logunit,*) " -- Program stopping ..."
       call LDT_endrun
    end select 
    deallocate( lat_line, lon_line )

! -------------------------------------------------------------------
!     AGGREGATING FINE-SCALE GRIDS TO COARSER LIS OUTPUT GRID
! -------------------------------------------------------------------
    mi = subpnc*subpnr
    mo = LDT_rc%lnc(n)*LDT_rc%lnr(n)
    if( mi .ne. mo .and. LDT_rc%lc_gridtransform(n) == "none" ) then
      write(LDT_logunit,*) "ERR: Spatial transform, 'none', is selected, but number of"
      write(LDT_logunit,*) "      input and output points do not match. Select other spatial"
      write(LDT_logunit,*) "      option (e.g., mode, neighbor, tile, etc.)."
      write(LDT_logunit,*) "Program stopping ..."
      call LDT_endrun
    endif

 !- No spatial transform applied:
    if( LDT_rc%lc_gridtransform(n) == "none" .and. &
        LDT_rc%lis_map_proj(n) == "latlon" ) then

       write(LDT_logunit,*) " No aggregation applied for parameter file ... "
       vegtype(:,:) = read_inputparm(:,:)

 !- Other transformations applied:
    else
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
       elseif( LDT_rc%lc_gridtransform(n) == "tile" ) then

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

end subroutine read_USGSNative_lc

