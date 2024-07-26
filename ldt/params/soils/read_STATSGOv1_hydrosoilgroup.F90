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
! !ROUTINE: read_STATSGOv1_hydrosoilgroup
! \label{read_STATSGOv1_hydrosoilgroup}
!
! !REVISION HISTORY:
!  03 Sept 2004: Sujay Kumar; Initial Specification
!  03  Oct 2013: KR Arsenault;  Expanded to read in native STATSGO v1 files
!
! !INTERFACE:
subroutine read_STATSGOv1_hydrosoilgroup( n, fgrd )

! !USES:
  use LDT_coreMod, only : LDT_rc
  use LDT_logMod,  only : LDT_logunit, LDT_getNextUnitNumber, &
                          LDT_releaseUnitNumber, LDT_endrun
  use LDT_fileIOMod
  use LDT_binaryIOlayer_module, only : LDT_read_bsqfile

  implicit none
! !ARGUMENTS: 
  integer,intent(in)    :: n
  real,   intent(inout) :: fgrd(LDT_rc%lnc(n),LDT_rc%lnr(n),5)
!
! !DESCRIPTION:
!  This subroutine retrieves STATSGO v1 hydrological soils group data and 
!  reprojects it to the latlon projection. 
!
!  The 30 sec soil data starts at: latitude 90. and longitude -180.
!
!    - the file is in direct access format
!    - each direct access record contains data in one latitude circle
!      beginning at -180 degree longitude, and end at +180 degree
!    - the data is arranged to start at northernmost latitude (north pole),
!      and end at south pole
!
!  The arguments are:
!  \begin{description}
!  \item[n]
!   index of the nest
!  \item[fgrd]
!   output field with the retrieved hydrological soils group
!  \end{description}
!EOP

   integer, parameter :: input_cols = 6936
   integer, parameter :: input_rows = 2984
   real,    parameter :: input_xres = 1.0/120.0
   real,    parameter :: input_yres = 1.0/120.0

   integer*4 :: hsg_types
 
   integer :: c, r, i, t 
   integer :: ftn
   integer :: length, nrec
   integer :: water_class
   logical :: file_exists
   integer :: mi                       ! Total number of input param grid fgrd points
   integer :: mo                       ! Total number of output LIS grid fgrd points
   real    :: param_gridDesc(20)       ! Input parameter grid desc fgrd
   integer, allocatable  :: n11(:)     ! Array that maps the location of each input grid
                                      !   point in the output grid. 
   real,    allocatable  :: gi(:)      ! Input parameter 1d grid
   logical*1,allocatable :: li(:)      ! Input logical mask (to match gi)
   real      :: go1(LDT_rc%lnc(n)*LDT_rc%lnr(n))            ! Output lis 1d grid
   logical*1 :: lo1(LDT_rc%lnc(n)*LDT_rc%lnr(n))            ! Output logical mask (to match go)
   real      :: temp(LDT_rc%lnc(n),LDT_rc%lnr(n))

   integer*2, allocatable :: soil_int(:,:,:)
   real,      allocatable :: input_hsg(:,:,:)

   character*50 :: hsg_proj
   character*50 :: hsg_gridtransform 

! ___________________________________________________________________
  
   water_class = 14    !--> ?? Water class for STATSGO v1 hydrological soils group class
   temp = float(water_class)
   fgrd = 0.

!- Set parameter grid fgrd inputs:
!  Below grid based on info from:
!   http://dbwww.essc.psu.edu/dbtop/.link/1995-0795/proj_geo.html
   hsg_proj = "latlon"
   param_gridDesc(1)  = 0.                ! Latlon
   param_gridDesc(2)  = input_cols
   param_gridDesc(3)  = input_rows
   param_gridDesc(4)  =  24.5330+(1/240)  ! LL lat
   param_gridDesc(5)  = -124.750+(1/240)  ! LL lon
   param_gridDesc(6)  = 128
   param_gridDesc(7)  =  49.400-(1/240)   ! UR lat
   param_gridDesc(8)  = -66.950-(1/240)   ! UR lon
   param_gridDesc(9)  = input_yres        ! dy: 0.0083333
   param_gridDesc(10) = input_xres        ! dx: 0.0083333
   param_gridDesc(20) = 64
!   param_gridDesc(20) = 255

   LDT_rc%hsg_proj = hsg_proj
   LDT_rc%hsg_gridDesc(n,:) = param_gridDesc(:)

   inquire(file=trim(LDT_rc%hsgfile(n)), exist=file_exists)
   if(.not.file_exists) then 
      write(LDT_logunit,*) "Hydrological soil group map ",&
                            trim(LDT_rc%hsgfile(n))," not found."
      write(LDT_logunit,*) "Program stopping ..."
      call LDT_endrun
   endif

   write(unit=LDT_logunit,fmt=*) "[INFO] Reading STATSGO v1 hydrological soil group file: ",&
        trim(LDT_rc%hsgfile(n))

 ! Check if Parameter domain is within LIS Run Domain:
   call LDT_checkDomainExtents(n, param_gridDesc(:))

!- Select gridtransform type by LIS and param grid resolutions:
   hsg_gridtransform = "none"

   if( LDT_rc%gridDesc(n,9)  == param_gridDesc(9) .and. &
       LDT_rc%gridDesc(n,10) == param_gridDesc(10) ) then
      hsg_gridtransform = "none"
   endif
   if( LDT_rc%gridDesc(n,9)  > param_gridDesc(9) .or. &
       LDT_rc%gridDesc(n,10) > param_gridDesc(10) ) then
      hsg_gridtransform = "mode"
   endif
!   if( LDT_rc%gridDesc(n,9)  < param_gridDesc(9) .or. &
!       LDT_rc%gridDesc(n,10) < param_gridDesc(10) ) then
!      hsg_gridtransform = "neighbor"
!   endif

   write(unit=LDT_logunit,fmt=*) "[INFO] HSG Grid transform is: ",&
         trim(hsg_gridtransform)

   hsg_types = 5    ! Number of HSG types 

   allocate( soil_int(input_cols, input_rows, hsg_types) )
   soil_int = water_class

!- Open the BSQ file:
   ftn = LDT_getNextUnitNumber()
   open( ftn, file=LDT_rc%hsgfile(n), &
         access="direct", form="formatted", recl=input_cols*1)
!  inbyte = 1

!- Read-in the BSQ file:
   call LDT_read_bsqfile( ftn, hsg_types, input_rows, input_cols, &
            1, 1, 1, 1, hsg_types, input_rows, input_cols, soil_int(:,:,:) )

   call LDT_releaseUnitNumber(ftn)

   allocate( input_hsg(input_cols, input_rows, hsg_types) )
   input_hsg = float(water_class)

!- Flip y-axis rows:
   i = 0
   do r = input_rows, 1, -1;  i = i + 1
      do c = 1, input_cols
         input_hsg(c,i,:) = float(soil_int(c,r,:))
      end do
   end do
   deallocate( soil_int )

   fgrd(:,:,:) = input_hsg(:,:,:)

   return

! -------------------------------------------------------------------
!     AGGREGATING FINE-SCALE GRIDS TO COARSER LIS OUTPUT GRID
! -------------------------------------------------------------------

#if 0
   mi = input_cols*input_rows
   mo = LDT_rc%lnc(n)*LDT_rc%lnr(n)
!   if( mi .ne. mo .and. LDT_rc%soils_gridtransform(n) == "none" ) then
   if( mi .ne. mo .and. hsg_gridtransform  == "none" ) then
      write(LDT_logunit,*) "[ERR] Spatial transform, 'none', is selected, but number of"
      write(LDT_logunit,*) "  input and output points do not match. Select other spatial"
      write(LDT_logunit,*) "  option (e.g., mode, tile, etc.)."
      write(LDT_logunit,*) "Program stopping ..."
      call LDT_endrun
   endif

   allocate( li(mi), gi(mi) )
   gi = float(LDT_rc%waterclass)
   li = .false.
   lo1 = .false.
!- Assign 2-D input soil text to 1-D for aggregation routines:
   i = 0
   do r = 1, input_rows
      do c = 1, input_cols;  i = i + 1
         gi(i) = input_hsg(c,r)
         if( gi(i) .ne. LDT_rc%udef ) li(i) = .true.
      enddo
   enddo

!- Apply the spatial transform option:
!   select case( LDT_rc%soils_gridtransform(n) )
   select case( hsg_gridtransform )

  !- (a) Single-layer selection:
     case( "none", "mode", "neighbor" )

     !- Transform parameter from original grid to LIS output grid:
!        call LDT_transform_paramgrid(n, LDT_rc%soils_gridtransform(n), &
        call LDT_transform_paramgrid(n,  hsg_gridtransform, &
                           param_gridDesc, mi, 1, gi, li, mo, go1, lo1 )

     !- Convert 1D count to 2D grid fgrds:
        i = 0
        do r = 1, LDT_rc%lnr(n)
           do c = 1, LDT_rc%lnc(n)
              i = i + 1
              temp(c,r) = go1(i)
           enddo
        enddo

       deallocate( gi, li )

   end select  ! End grid cnt aggregation method
   deallocate( input_hsg )

!-  2-D Array space:
!   if( LDT_rc%soils_gridtransform(n) == "none" .or. &  ! Non-tiled surfaces
!       LDT_rc%soils_gridtransform(n) == "mode" .or. &
!       LDT_rc%soils_gridtransform(n) == "average" ) then
   select case( hsg_gridtransform )

    case( "none" )

   !- Assign hydrological soils group types of less than 0 to an actual hsg value:
      do r = 1, LDT_rc%lnr(n)
         do c = 1, LDT_rc%lnc(n)
            if( nint(temp(c,r)) .le. 0 ) then
              temp(c,r) = 12   ! STATSGO -- Clay 
            endif
!            if( (nint(temp(c,r)) .ne. water_class  ) .and. &
!                (nint(temp(c,r)) .ne. LDT_rc%udef) ) then
!               gridcnt(c,r,NINT(temp(c,r))) = 1.0
!            endif
         enddo
      enddo

   end select
!   end if

#endif

   write(LDT_logunit,fmt=*) "[INFO] Done reading STATSGO v1 hydrological soil group file"

 end subroutine read_STATSGOv1_hydrosoilgroup

