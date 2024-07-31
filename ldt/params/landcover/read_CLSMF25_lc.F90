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
! !ROUTINE: read_CLSMF25_lc
!  \label{read_CLSMF25_lc}
!
! !REVISION HISTORY:
!  03 Sept 2004: Sujay Kumar; Initial Specification
!  23  May 2012: KR Arsenault; Implemented new features to read different
!                 resolutions and generate landmask from landcover
!
! !INTERFACE:
subroutine read_CLSMF25_lc(n, num_types, fgrd, maskarray)

! !USES:
  use LDT_coreMod, only : LDT_rc, LDT_domain
  use LDT_logMod,  only : LDT_logunit, LDT_getNextUnitNumber, &
                          LDT_releaseUnitNumber, LDT_verify, LDT_endrun
  use map_utils

  implicit none

! !ARGUMENTS: 
  integer, intent(in) :: n
  integer, intent(in) :: num_types
  real, intent(inout) :: fgrd(LDT_rc%lnc(n),LDT_rc%lnr(n),LDT_rc%nt)
  real, intent(inout) :: maskarray(LDT_rc%lnc(n),LDT_rc%lnr(n))
!
! !DESCRIPTION:
!  This subroutine reads the AVHRR landcover data and returns the 
!  distribution of vegetation in each grid cell, in a lat/lon
!  projection.  Also, the landmask is either generated and/or 
!  read in this routine.
!
! Mosaic/CLSM F2.5 Land cover classes:
!  1 - Broadleaf evergreen trees
!  2 - Broadleaf deciduous trees
!  3 - Needleleaf trees
!  4 - Grassland
!  5 - Broadleaf shrubs
!  6 - Dwarf trees
!
! ** Only 6 types in the Catchment Fortuna-2.5 and later.
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
   integer :: i, t, c, r, line, gc, gr
   integer :: dummy_int
   integer :: glpnr, glpnc
   real    :: param_grid(20)
   real    :: rlat(LDT_rc%lnc(n),LDT_rc%lnr(n))
   real    :: rlon(LDT_rc%lnc(n),LDT_rc%lnr(n))
   integer, allocatable :: tmpint(:), tmptileid(:)
   real,    allocatable :: vegtype(:,:)
!__________________________________________________________________

   maskarray(:,:) = 0.0
   fgrd(:,:,:)    = 0.0

!- Determine global/complete parameter domain number of points:
   param_grid(:) = LDT_rc%mask_gridDesc(n,:)
   glpnr = nint((param_grid(7)-param_grid(4))/param_grid(10)) + 1
   glpnc = nint((param_grid(8)-param_grid(5))/param_grid(9)) + 1

!- Assign additional land cover types, including generic water points: 
   select case ( LDT_rc%lc_type(n) )
    case ( "UMD" ) 
      LDT_rc%bareclass    = 12
      LDT_rc%urbanclass   = 13
      LDT_rc%waterclass   = 14
      LDT_rc%snowclass    = 0
      LDT_rc%wetlandclass = 0
      LDT_rc%glacierclass = 0
    case ( "MOSAIC" )
      LDT_rc%waterclass   = 7
      LDT_rc%bareclass    = 5
      LDT_rc%urbanclass   = 5
      LDT_rc%glacierclass = 5
      LDT_rc%wetlandclass = 0
      LDT_rc%snowclass    = 0
    case default ! non-supported options
      write(LDT_logunit,*) "The land classification: ",trim(LDT_rc%lc_type(n)),&
                           " does not exist for CLSM F2.5 source."
      write(LDT_logunit,*) " -- Please select either: UMD or MOSAIC "
      write(LDT_logunit,*) "Program stopping ..."
      call LDT_endrun
   end select

!- Double-check landcover spatial transform option:
   if( LDT_rc%lc_gridtransform(n)=="tile" .or. LDT_rc%lc_gridtransform(n)=="mode")then
      write(LDT_logunit,*) " (in read_CLSMF25_lc) :: The 'tile' or 'mode' spatial transform option"
      write(LDT_logunit,*) "          has been selected, but these options are not"
      write(LDT_logunit,*) "          currently supported for the Catchment F2.5 LSM. "
      write(LDT_logunit,*) "          Please select option 'none' and run LDT again. "
      write(LDT_logunit,*) " Program stopping ..."
      call LDT_endrun
   end if

!- Check if land cover file exists:
   inquire( file=trim(LDT_rc%vfile(n)), exist=file_exists ) 
   if(.not. file_exists) then 
      write(LDT_logunit,*) "Landcover map: ",trim(LDT_rc%vfile(n))," does not exist."
      write(LDT_logunit,*) "Program stopping ..."
      call LDT_endrun
   endif

! -------------------------------------------------------------------
!    READ IN LANDMASK FOR CLSM ...
! -------------------------------------------------------------------

!- "READ-IN" land mask file, if user-specified:
   if( LDT_rc%mask_type(n) == "readin" ) then

     call read_clsm_maskfile(n, maskarray)

!- "CREATE" land mask and surface type fields (user-specified):
   else
      write(LDT_logunit,*) " (in read_CLSMF25_lc) :: The CLSM F2.5 parameters currently"
      write(LDT_logunit,*) "          require a defined land/water mask, so a mask cannot"
      write(LDT_logunit,*) "          be 'created' or any other option applied at this time. "
      write(LDT_logunit,*) "          Please select option 'readin' and run LDT again. "
      write(LDT_logunit,*) " Program stopping ..."
      call LDT_endrun
   end if

! -------------------------------------------------------------------
!    READ IN LAND COVER PARAMETER FIELDS (NON-TILED OPTION for now)
! -------------------------------------------------------------------

!- Open LDT land cover file:
   write(unit=LDT_logunit,fmt=*) "[INFO] Reading landcover file:", trim(LDT_rc%vfile(n))
   ftn = LDT_getNextUnitNumber()
   open(ftn, file=LDT_rc%vfile(n), form='formatted', status='old',iostat=ios1)

   allocate( tmpint(LDT_rc%nmaskpts(n)), tmptileid(LDT_rc%nmaskpts(n)) )
   tmptileid = 0
   tmpint = LDT_rc%waterclass

!- Read CLSM F2.5 vegclass (ascii) array:
   do t = 1, LDT_rc%nmaskpts(n)
      read(ftn,*)  tmptileid(t), dummy_int, tmpint(t)
   end do
 
!   allocate( vegtype(LDT_rc%lnc(n),LDT_rc%lnr(n)), stat=ierr )  ! Subsetted
   allocate( vegtype(glpnc,glpnr), stat=ierr )                   ! Complete domain
   call LDT_verify(ierr,'Error allocating vegtype')
   vegtype = float(LDT_rc%waterclass)

   t = 0
! - For future subsetted domains:
!   do r = 1, LDT_rc%lnr(n)
!      do c = 1, LDT_rc%lnc(n)
!         if( maskarray(c,r) > 0. ) then
! - For now - complete domains:
   do r = 1, glpnr
      do c = 1, glpnc

         if( LDT_rc%global_mask(c,r) > 0. ) then
           t = t + 1
           if( tmpint(t) > 0 ) then
              vegtype(c,r) = float(tmpint(t))
           elseif( tmpint(t) == 0 ) then
              vegtype(c,r) = 5.   ! Assign landcover for glacier points
           endif
         endif

      end do
   enddo
!   write(*,*) count( mask = vegtype(:,:).ne. 7)

! -------------------------------------------------------------------
!    CREATE OR READ-IN (OR IMPOSE) LAND MASK FILE AND CREATE
!    SURFACE MAP
! -------------------------------------------------------------------

!!- Final fgrd output fields:
 !- Build in subsetting domain:
   do r = 1, LDT_rc%lnr(n)
      do c = 1, LDT_rc%lnc(n)

         call ij_to_latlon(LDT_domain(n)%ldtproj,float(c),float(r),&
                           rlat(c,r),rlon(c,r))
         gr = nint((rlat(c,r)-param_grid(4))/param_grid(10))+1
         gc = nint((rlon(c,r)-param_grid(5))/param_grid(9))+1

         if( vegtype(gc,gr) > 0. .and. vegtype(gc,gr).ne.LDT_rc%waterclass ) then 
            fgrd(c,r,int(vegtype(gc,gr))) = 1.0
         else
            fgrd(c,r,LDT_rc%waterclass) = 1.0
         endif
      enddo
   enddo
   
   deallocate( tmpint, tmptileid )
   deallocate( vegtype )

   call LDT_releaseUnitNumber(ftn)

end subroutine read_CLSMF25_lc


!BOP
!
! !ROUTINE: read_clsm_maskfile
!  \label{read_clsm_maskfile}
!
! !REVISION HISTORY:
!  03 Sept 2004: Sujay Kumar; Initial Specification
!  01 June 2012: KR Arsenault; Restructured to simply read in a mask file
!
! !INTERFACE:
 subroutine read_clsm_maskfile( n, localmask )

! !USES:
  use LDT_coreMod, only : LDT_rc, LDT_localPet
  use LDT_logMod,  only : LDT_logunit, LDT_getNextUnitNumber, &
                          LDT_releaseUnitNumber, LDT_endrun
  use LDT_gridmappingMod

  implicit none

! !ARGUMENTS: 
  integer, intent(in)  :: n
  real,    intent(out) :: localmask(LDT_rc%lnc(n),LDT_rc%lnr(n))

! !DESCRIPTION:
!  This subroutine reads the landmask data and returns the 
!   mask and surface type arrays.
!
!  The arguments are:
!  \begin{description}
!   \item[n]
!    index of nest
!   \item[localmask]
!    landmask for the region of interest
!   \end{description}
!
!EOP      
  integer :: ftn, ios1
  logical :: file_exists
  integer :: c, r, t, line
  integer :: glpnc, glpnr             ! Parameter (global) total columns and rows
  integer :: subpnc, subpnr           ! Parameter subsetted columns and rows
  real    :: subparam_gridDesc(20)    ! Input parameter grid desc array

  integer, allocatable :: lat_line(:,:)
  integer, allocatable :: lon_line(:,:)
  real,    allocatable :: read_inputparm(:,:)  ! Read input parameter
!_________________________________________________________________________________

   LDT_rc%nmaskpts = 0.
   localmask = 0.

!- Check for and open landmask file:
   inquire(file=trim(LDT_rc%mfile(n)), exist=file_exists)
   if( file_exists ) then 
      write(LDT_logunit,*)'[INFO] Reading CLSM mask file:',trim(LDT_rc%mfile(n)), & 
                          ' (',LDT_localPet,')'

! -------------------------------------------------------------------
!    PREPARE SUBSETTED PARAMETER GRID FOR READING IN NEEDED DATA
! -------------------------------------------------------------------
 !- Map Parameter Grid Info to LIS Target Grid/Projection Info -- 
    subparam_gridDesc = 0.

    call LDT_RunDomainPts( n, LDT_rc%mask_proj, LDT_rc%mask_gridDesc(n,:), &
                  glpnc, glpnr, subpnc, subpnr,  &
                  subparam_gridDesc, lat_line, lon_line )

    allocate( read_inputparm(subpnc, subpnr) )
    read_inputparm = 0.

 ! -------------------------------------------------------------------

 ! -- Open land/water mask file:
      ftn = LDT_getNextUnitNumber()
      open(ftn, file=trim(LDT_rc%mfile(n)),form='unformatted', recl=4, &
           access='direct', iostat=ios1)

 ! == (1) READ IN GLOBAL/ENTIRE MASK PARAMETER DATA: == 
 !     (Currently important for reading in other CLSM F2.5 parameters)

    ! Global mask field for current CLSM F2.5 needs:
      allocate(LDT_rc%global_mask(glpnc,glpnr))  !  Temporary ...
      LDT_rc%global_mask = LDT_rc%udef
      line = 0
      do r = 1, glpnr 
         do c = 1, glpnc
            line = line + 1
            read(ftn,rec=line)  LDT_rc%global_mask(c,r)
         enddo
      enddo
 !  ( to be removed later after full CLSM-F2.5 preprocesser
 !    is implemented into LDT)
 ! =======================================================

 ! == (2) READ IN LANDMASK DATA: == 
      line = 0
      do r = 1, subpnr
         do c = 1, subpnc
            line = (lat_line(c,r)-1)*glpnc + lon_line(c,r)
            read(ftn,rec=line) read_inputparm(c,r)
         enddo
      enddo
      localmask = read_inputparm

 ! == (3) INCLUDE WATER POINTS, IF SELECTED: == 

      if( LDT_rc%inc_water_pts ) then
         write(*,*) " FOR CLSM F2.5, PLEASE DO NOT INCLUDE WATER PTS "
         write(*,*) "  AT THIS TIME .... stopping."
         call LDT_endrun

         do r = 1, glpnr
            do c = 1, glpnc
               if( LDT_rc%global_mask(c,r) == 0 .or. &
                   LDT_rc%global_mask(c,r) == LDT_rc%waterclass ) then
                 LDT_rc%global_mask(c,r) = 1 
               endif
            enddo
         enddo
         do r = 1, LDT_rc%lnr(n)
            do c = 1, LDT_rc%lnc(n)
               if( localmask(c,r) == 0 .or. &
                   localmask(c,r) == LDT_rc%waterclass ) then
                 localmask(c,r) = 1 
               endif
            end do
         end do
      end if

   !- Generate total number of accounted mask points:
! - Eventually will handle subsetted mask:
!      do r=1,LDT_rc%lnr(n)
!         do c=1,LDT_rc%lnc(n)
!             if( localmask(c,r) >= 1 ) then
! - Set for global mask for now:
      do r=1,glpnr
         do c=1,glpnc
            if( LDT_rc%global_mask(c,r) >= 1 ) then
                LDT_rc%nmaskpts(n) = LDT_rc%nmaskpts(n) + 1
            endif
         end do 
      end do

      call LDT_releaseUnitNumber(ftn)

   else
      write(LDT_logunit,*) "Landmask map: ",trim(LDT_rc%mfile(n)), " does not exist"
      write(LDT_logunit,*) "program stopping ..."
      call LDT_endrun
   endif

end subroutine read_clsm_maskfile
