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
! !ROUTINE: read_AVHRR_GFS_lc
!  \label{read_AVHRR_GFS_lc}
!
! !REVISION HISTORY:
!  03 Sept 2004: Sujay Kumar; Initial Specification
!  23  May 2012: KR Arsenault; Implemented new features to read different
!                 resolutions and generate landmask from landcover
!
! !INTERFACE:
subroutine read_AVHRR_GFS_lc(n, num_types, fgrd, maskarray)

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
!  This subroutine reads the GFS landcover data (from AVHRR) and returns the 
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
   integer :: t, c, r
   real    :: isum
   integer :: glpnc, glpnr             ! Parameter (global) total columns and rows
   integer :: subpnc, subpnr           ! Parameter subsetted columns and rows
   real    :: subparam_gridDesc(20)    ! Input parameter grid desc array
   integer, allocatable  :: lat_line(:,:), lon_line(:,:)
   real,    allocatable  :: dummy_lat(:,:), dummy_lon(:,:), whole_veg_layer(:,:)
   real,    allocatable  :: vegtype(:,:)
   real,    allocatable  :: read_vegcnt(:,:,:)   ! Read input parameter
!__________________________________________________________________

!- Assign additional land cover types, including generic water points: 
   select case (LDT_rc%lc_type(n))
    case ("UMD") 
      LDT_rc%bareclass    = 12
      LDT_rc%urbanclass   = 13
      LDT_rc%waterclass   = 14
      LDT_rc%snowclass    = 0
      LDT_rc%wetlandclass = 0
      LDT_rc%glacierclass = 0
    case default ! non-supported options
      write(LDT_logunit,*) "The land classification: ",trim(LDT_rc%lc_type(n)),&
                           " does not exist for AVHRR_GFS source."
      write(LDT_logunit,*) " -- Please select UMD"
      write(LDT_logunit,*) "Program stopping ..."
      call LDT_endrun
   end select

!- Check if land cover file exists:
   inquire(file=trim(LDT_rc%vfile(n)), exist=file_exists) 
   if(.not. file_exists) then 
      write(LDT_logunit,*) "Landcover map: ",trim(LDT_rc%vfile(n))," does not exist."
      write(LDT_logunit,*) "Program stopping ..."
      call LDT_endrun
   endif
   write(unit=LDT_logunit,fmt=*)"[INFO] Reading landcover file:", &
       trim(LDT_rc%vfile(n))

!- Open LDT land cover file:
   ftn = LDT_getNextUnitNumber()
   open(ftn, file=LDT_rc%vfile(n),status='old',form='unformatted',&
       access='sequential',iostat=ios1)
! -------------------------------------------------------------------
   maskarray = 0.0
   fgrd      = 0.0
! -------------------------------------------------------------------
!- Map Parameter Grid Info to LIS Target Grid/Projection Info -- 
   subparam_gridDesc = 0.
   call LDT_RunDomainPts( n, LDT_rc%lc_proj, LDT_rc%lc_gridDesc(n,:), &
                          glpnc, glpnr, subpnc, subpnr,  &
                          subparam_gridDesc, lat_line, lon_line )
! -------------------------------------------------------------------
   write(unit=LDT_logunit,fmt=*) "[INFO] Reading a 3-D landcover map ..."

   if(LDT_rc%lc_gridtransform(n).ne."none".and. &
      LDT_rc%lc_gridtransform(n).ne.'tile') then
        write(LDT_logunit,*) "[ERR] Landcover file has gaussian grid coordinate system and"
        write(LDT_logunit,*) "  therefore cannot be translated to different output grid and "
        write(LDT_logunit,*) "  projection at this time."
        write(LDT_logunit,*) "  Must select the 'none' or 'tile' spatial transform option "
        write(LDT_logunit,*) "  (in this case both options do the same thing)."
        write(LDT_logunit,*) "Program stopping ..."
        call LDT_endrun
   endif

   if(LDT_rc%lis_map_proj(n).ne."gaussian") then
     write(LDT_logunit,*) "[ERR] Reading in GFS UMD 3D tiled veg file"   
     write(LDT_logunit,*) "[ERR] only gaussian LIS grids are supported"
     write(LDT_logunit,*) "     projection: ",trim(LDT_rc%lis_map_proj(n))
     write(LDT_logunit,*) " Stopping run ..."
     call LDT_endrun
   endif
! -------------------------------------------------------------------
!    READ IN LAND COVER PARAMETER FIELDS (NON-TILED/TILED OPTIONS)
! -------------------------------------------------------------------
   select case( LDT_rc%lc_gridtransform(n) ) 
     case( "none", "tile" ) 
       allocate( read_vegcnt(subpnc,subpnr,LDT_rc%nt) )
       read_vegcnt = 0.

       allocate(dummy_lat(glpnc,glpnr))
       allocate(dummy_lon(glpnc,glpnr))
       allocate(whole_veg_layer(glpnc,glpnr))
       read(ftn) dummy_lat
       read(ftn) dummy_lon
       do t = 1,LDT_rc%nt
         if ( t == LDT_rc%nt .and. glpnr == 190 ) then
            ! NR = 190 corresponds to 95 latitude circles and the T126 domain.
            ! This domain contains 13 veg types not 14.
            exit
         endif
         read(ftn) whole_veg_layer
         do r = 1,subpnr
           do c = 1,subpnc
              read_vegcnt(c,r,t) = whole_veg_layer(lon_line(c,r),lat_line(c,r))
           enddo ! columns
         enddo ! rows
       enddo ! tilespace
       deallocate(dummy_lat)
       deallocate(dummy_lon)
       deallocate(whole_veg_layer)
         
!       ! deal with water class
!       do r = 1,subpnr
!         do c = 1,subpnc
!           isum = sum(read_vegcnt(c,r,1:LDT_rc%nt))
!           read_vegcnt(c,r,LDT_rc%waterclass) = 1-isum
!         enddo
!      enddo
  endselect
  fgrd = read_vegcnt
  deallocate(read_vegcnt)
   
  ! create vegtype array necessary for landmask checks
  allocate(vegtype(LDT_rc%lnc(n),LDT_rc%lnr(n)))
  vegtype = -999
  do r = 1,subpnr
    do c = 1,subpnc
      do t = 1,LDT_rc%nt
        if (fgrd(c,r,t).gt.0.and.t.ne.LDT_rc%waterclass) vegtype(c,r) = t
        if (fgrd(c,r,t).gt.0.and.t.eq.LDT_rc%waterclass) vegtype(c,r) = 0
      enddo 
    enddo
  enddo
 
! -------------------------------------------------------------------
!    CREATE OR READ-IN (OR IMPOSE) LAND MASK FILE AND CREATE
!    SURFACE MAP
! -------------------------------------------------------------------
   if( LDT_rc%mask_type(n) == "readin" ) then
     call read_gfs_maskfile(n,vegtype,maskarray)
   elseif( LDT_rc%mask_type(n) == "create" ) then
     write(LDT_logunit,*) "[ERR] Create Landmask option not supported for Gaussian projection at this time" 
     write(LDT_logunit,*) "Program stopping ..."
     call LDT_endrun
   endif

   call LDT_releaseUnitNumber(ftn)
end subroutine read_AVHRR_GFS_lc


subroutine read_gfs_maskfile(n,vegtype,array)
  use LDT_coreMod,     only : LDT_rc, LDT_domain
  use LDT_logMod,      only : LDT_logunit, LDT_getNextUnitNumber, &
          LDT_releaseUnitNumber, LDT_endrun
  use LDT_gridmappingMod

  implicit none
  integer, intent(in) :: n
  real, intent(in) :: vegtype(LDT_rc%lnc(n),LDT_rc%lnr(n))
  real, intent(out) :: array(LDT_rc%lnc(n),LDT_rc%lnr(n))

  integer :: ftn
  integer :: c,r,t
  logical :: file_exists
  integer :: subpnc, subpnr, glpnc, glpnr
  real    :: subparam_gridDesc(20)       
  integer, allocatable :: lat_line(:,:), lon_line(:,:)
  real,    allocatable :: read_mask(:,:), dummy_lat(:,:), dummy_lon(:,:)

! check
  inquire(file=trim(LDT_rc%mFile(n)), exist=file_exists)
  if(.not.file_exists) then
     write(LDT_logunit,*) "Landmask map ",trim(LDT_rc%mFile(n))," not found."
     write(LDT_logunit,*) "Program stopping ..."
     call LDT_endrun
  endif

! open file
  ftn = LDT_getNextUnitNumber()
  open(ftn, file=trim(LDT_rc%mFile(n)), access='sequential',status='old', &
       form="unformatted", recl=4)

! segregate subdomain
  subparam_gridDesc = 0.
  call LDT_RunDomainPts(n,LDT_rc%mask_proj,LDT_rc%mask_gridDesc(n,:), &
                         glpnc,glpnr,subpnc,subpnr,  &
                         subparam_gridDesc,lat_line,lon_line )

! read source
  allocate(dummy_lat(glpnc,glpnr))
  allocate(dummy_lon(glpnc,glpnr))
  allocate(read_mask(glpnc,glpnr))

  read(ftn) dummy_lat
  read(ftn) dummy_lon
  read(ftn) read_mask

  close(ftn)
  deallocate(dummy_lat)
  deallocate(dummy_lon)

! subparam
  select case (LDT_rc%mask_gridtransform(n))
    case("none")
      do r = 1,subpnr
        do c = 1,subpnc
           array(c,r) = read_mask(lon_line(c,r),lat_line(c,r))
        enddo ! columns
      enddo ! rows
      deallocate(read_mask)
    case default
       write(LDT_logunit,*) "[ERR] The spatial transform, ",&
                            trim(LDT_rc%mask_gridtransform(n)),&
                            ", for Landmask is not available at this time ... "
       write(LDT_logunit,*) "Program stopping ..."
       call LDT_endrun
   end select

! check against landcover
  t = 0
  do r = 1,subpnr
    do c = 1,subpnc
      if (array(c,r).eq.1.and.vegtype(c,r).le.0) then
        t = t+1
        write(LDT_logunit,*) "[ERR] Landmask and veg classes do not agree, ",&
               lon_line(c,r),lat_line(c,r),array(c,r),vegtype(c,r),t
!       write(LDT_logunit,*) "Program stopping ..."
!       call LDT_endrun
     endif 
   enddo
 enddo
  call LDT_releaseUnitNumber(ftn)
end subroutine read_gfs_maskfile



