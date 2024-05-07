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
! !ROUTINE: read_CLM45_lc
!  \label{read_CLM45_lc}
!
! !REVISION HISTORY:
!  21  Dec 2016: H Beaudoing; Initial Specification
!
! !INTERFACE:
subroutine read_CLM45_lc(n, num_types, fgrd, maskarray)

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
!  This subroutine reads the CLM-4.5 surface data and returns the 
!  distribution of landcover in each grid cell, in a lat/lon
!  projection.  Also, the landmask is either generated and/or 
!  read in this routine.
!
!max_pft_per_gcell = numpft+1 + 3 + maxpatch_urb*numurbl + maxpatch_glcmec
! numpft  maxpatch_urb numurbl  maxpatch_glcmec  numcft max_pft_per_gcell
! ---------------------------------------------------------------------------
!    16         5         3         0         2        35
!YDT  intercepting pft data and save to (overwrite) LDT file, run with single CPU
! clm pft 0-16: map to landcover and surfacetype(:, :, z), z=1-17
!     pft 61-65 (urban):  z: 18-22, 23-27, 28-32; (max 3 urban landunits )
!     deep lake:  3 (itype_lun)  z: 33
!     wetland:  5 (itype_lun)  z: 34
!     glacier: 2 (itype_lun)   z: 35
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
   logical :: file_exists
   integer :: i, t, c, r, line, gc, gr
   integer :: dummy_int
   integer :: glpnr, glpnc
   real    :: param_grid(20)
   real    :: rlat(LDT_rc%lnc(n),LDT_rc%lnr(n))
   real    :: rlon(LDT_rc%lnc(n),LDT_rc%lnr(n))
   real*8, allocatable :: tmpint(:,:,:)
!__________________________________________________________________

   maskarray(:,:) = 0.0
   fgrd(:,:,:)    = 0.0

!- Determine global/complete parameter domain number of points:
   param_grid(:) = LDT_rc%mask_gridDesc(n,:)
   glpnr = nint((param_grid(7)-param_grid(4))/param_grid(10)) + 1
   glpnc = nint((param_grid(8)-param_grid(5))/param_grid(9)) + 1

!- Assign additional land cover types, including generic water points: 
   select case ( LDT_rc%lc_type(n) )
    case ( "CLM45" ) 
      LDT_rc%bareclass    = 1 
      LDT_rc%urbanclass   = 18 !18-22, 23-27,28-32
      LDT_rc%waterclass   = 36
      LDT_rc%lakeclass    = 33
      LDT_rc%wetlandclass = 34
      LDT_rc%glacierclass = 35
    case default ! non-supported options
      write(LDT_logunit,*) "The land classification: ",trim(LDT_rc%lc_type(n)),&
                           " does not exist for CLM-4.5 source."
      write(LDT_logunit,*) " -- Please select : CLM45 "
      write(LDT_logunit,*) "Program stopping ..."
      call LDT_endrun
   end select

!- Double-check landcover spatial transform option:
   if( LDT_rc%lc_gridtransform(n)=="tile" .or. LDT_rc%lc_gridtransform(n)=="mode")then
      write(LDT_logunit,*) " (in read_CLM45_lc) :: The 'tile' or 'mode' spatial transform option"
      write(LDT_logunit,*) "          has been selected, but these options are not"
      write(LDT_logunit,*) "          currently supported for the CLM-4.5 LSM. "
      write(LDT_logunit,*) "          Please select option 'none' and run LDT again. "
      write(LDT_logunit,*) " Program stopping ..."
      call LDT_endrun
   end if

!- Check if land cover file exists:
   inquire( file=trim(LDT_rc%vfile(n)), exist=file_exists ) 
   if(.not. file_exists) then 
      write(LDT_logunit,*) "Landcover map: ",trim(LDT_rc%vfile(n)),&
                           " does not exist."
      write(LDT_logunit,*) "Program stopping ..."
      call LDT_endrun
   endif

! -------------------------------------------------------------------
!    READ IN LANDMASK FOR CLM-4.5 ...
! -------------------------------------------------------------------

!- "READ-IN" land mask file, if user-specified:
   if( LDT_rc%mask_type(n) == "readin" ) then

     call read_clm45_maskfile(n, glpnc, glpnr, maskarray)

!- "CREATE" land mask and surface type fields (user-specified):
   else
      write(LDT_logunit,*) " (in read_CLM45_lc) ::"
      write(LDT_logunit,*) "  The CLM-4.5 parameters currently require"
      write(LDT_logunit,*) "      a defined land/water mask, so a mask"
      write(LDT_logunit,*) "      cannot be 'created' or any other"
      write(LDT_logunit,*) "      option applied at this time. "
      write(LDT_logunit,*) "  Please select option 'readin'"
      write(LDT_logunit,*) "         and run LDT again. "
      write(LDT_logunit,*) " Program stopping ..."
      call LDT_endrun
   end if
   write(LDT_logunit,*) "maskarray values.",maxval(maskarray(:,:)),minval(maskarray(:,:))

! -------------------------------------------------------------------
!    GET LAND COVER PARAMETER FIELDS 
!    This is just for place holder.  Actual Landcover and Surfacetype
!    will be overwritten based on CLM-45 surface data in 
!    params/CLM45/CLM45_parmsMod.F90 routine
! -------------------------------------------------------------------
   allocate( tmpint(LDT_rc%lnc(n),LDT_rc%lnr(n),LDT_rc%nt) )
   tmpint = 0.0

   do r = 1, LDT_rc%lnr(n)
      do c = 1, LDT_rc%lnc(n)
         if( maskarray(c,r) > 0. ) then  ! using CLM-4.5 mask
           do t = 1, LDT_rc%nt    ! including 36=water
            fgrd(c,r,t) = real(tmpint(c,r,t))
           enddo
         else
            fgrd(c,r,LDT_rc%waterclass) = 1.0
         endif

      end do
   enddo
   
   deallocate( tmpint )

end subroutine read_CLM45_lc


!BOP
!
! !ROUTINE: read_clm45_maskfile
!  \label{read_clm45_maskfile}
!
! !INTERFACE:
 subroutine read_clm45_maskfile( n, glpnc, glpnr, localmask )

! !USES:
  use LDT_coreMod, only : LDT_rc, LDT_localPet
  use LDT_logMod,  only : LDT_logunit, LDT_getNextUnitNumber, &
                          LDT_releaseUnitNumber, LDT_endrun
  use LDT_gridmappingMod
  use CLM45_parmsMod

  implicit none

! !ARGUMENTS: 
  integer, intent(in)  :: n
  integer, intent(in)  :: glpnc, glpnr             ! Parameter (global) total columns and rows
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
  logical :: file_exists
  integer :: c, r, t, line
  integer :: glpnc2, glpnr2           ! global columns and rows dummy
  integer :: subpnc, subpnr           ! Parameter subsetted columns and rows
  real    :: subparam_gridDesc(20)    ! Input parameter grid desc array
  real    :: shifted_gridDesc(20)     ! Input parameter grid desc 180W-180E lon

  integer, allocatable :: lat_line(:,:)
  integer, allocatable :: lon_line(:,:)
  real*8, allocatable  :: dvalue(:,:,:,:)
!_________________________________________________________________________________

   LDT_rc%nmaskpts = 0.
   localmask = 0.

!- Check for and open landmask file:
   inquire(file=trim(LDT_rc%mfile(n)), exist=file_exists)
   if( file_exists ) then 
      write(LDT_logunit,*)'[INFO] Reading CLM45 mask file:',trim(LDT_rc%mfile(n)), & 
                          ' (',LDT_localPet,')'


 ! -------------------------------------------------------------------
 ! == (1) READ IN GLOBAL/ENTIRE MASK PARAMETER DATA: == 
 ! -------------------------------------------------------------------

    ! Global mask field for current CLM-4.5 needs:
      allocate(LDT_rc%global_mask(glpnc,glpnr))  !  Temporary ...
      allocate(dvalue(glpnc,glpnr,1,1))  !  Temporary ...
      LDT_rc%global_mask = LDT_rc%udef
      if( LDT_rc%lsm == "CLM.4.5") then
        call read_netcdf_4d_global(n,trim(LDT_rc%mfile(n)),"mask",2, &
             "nj","ni","none","none",glpnc,glpnr,1,1,&
             dvalue(:,:,1,1),shifted_gridDesc)
        write(LDT_logunit,*) "Done reading CLM-4.5 mask values.",maxval(dvalue(:,:,1,1)),minval(dvalue(:,:,1,1))
        write(LDT_logunit,*) "shifted_gridDesc= ",shifted_gridDesc
        do r = 1, glpnr
           do c = 1, glpnc
              LDT_rc%global_mask(c,r) = real(dvalue(c,r,1,1))
           enddo
        enddo
      endif
      deallocate(dvalue)
      write(LDT_logunit,*) "global_mask mask values.",maxval(LDT_rc%global_mask(:,:)),minval(LDT_rc%global_mask(:,:))

 ! -------------------------------------------------------------------
 !    PREPARE SUBSETTED PARAMETER GRID FOR READING IN NEEDED DATA
 ! -------------------------------------------------------------------
 !- Map Parameter Grid Info to LIS Target Grid/Projection Info -- 
    subparam_gridDesc = 0.
    glpnc2 = glpnc
    glpnr2 = glpnr
    call LDT_RunDomainPts( n, LDT_rc%mask_proj, shifted_gridDesc(:), &
                  glpnc2, glpnr2, subpnc, subpnr,  &
                  subparam_gridDesc, lat_line, lon_line )

 !- Subset mask field
    if ( glpnc2 == subpnc .and. glpnr2 == subpnr ) then
     localmask = LDT_rc%global_mask  ! NO regional mask option
    else
     do r = 1, subpnr
      do c = 1, subpnc
       localmask(c,r) = LDT_rc%global_mask(lon_line(c,r),lat_line(c, r))
!       print*,'here',c,r,lon_line(c,r),lat_line(c, r)
      enddo   !c
     enddo   !r
    endif   ! subset
    deallocate( lat_line, lon_line )

 ! == (2) Downscale/Aggregate mask points, if option turned on ==
 ! == (3) INCLUDE WATER POINTS, IF SELECTED: == 

      if( LDT_rc%inc_water_pts ) then
         write(*,*) " FOR CLM-4.5, PLEASE DO NOT INCLUDE WATER PTS"
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
      do r=1,LDT_rc%lnr(n)
         do c=1,LDT_rc%lnc(n)
            if( localmask(c,r) >= 1 ) then
                LDT_rc%nmaskpts(n) = LDT_rc%nmaskpts(n) + 1
            endif
         end do 
      end do
      write(LDT_logunit,*) "nmaskpts: ",LDT_rc%nmaskpts(n)

   else
      write(LDT_logunit,*) "Landmask map: ",trim(LDT_rc%mfile(n)),&
                           " does not exist"
      write(LDT_logunit,*) "program stopping ..."
      call LDT_endrun
   endif

end subroutine read_clm45_maskfile
