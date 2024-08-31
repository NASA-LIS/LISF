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
! !ROUTINE: read_SACHTET356_lc
!  \label{read_SACHTET356_lc}
!
! !REVISION HISTORY:
!  03 Sept 2004: Sujay Kumar; Initial Specification
!  23  May 2012: KR Arsenault; Implemented new features to read different
!                 resolutions and generate landmask from landcover
!
! !INTERFACE:
subroutine read_SACHTET356_lc(n, num_types, fgrd, maskarray)

! !USES:
  use LDT_coreMod,  only : LDT_rc
  use LDT_logMod,   only : LDT_logunit, LDT_getNextUnitNumber, &
          LDT_releaseUnitNumber, LDT_verify, LDT_endrun
  use LDT_gridmappingMod    
  use LDT_paramTileInputMod, only: param_index_fgrdcalc
  use LDT_xmrg_reader

  implicit none

! !ARGUMENTS: 
  integer, intent(in) :: n
  integer, intent(in) :: num_types
  real, intent(inout) :: fgrd(LDT_rc%lnc(n),LDT_rc%lnr(n),LDT_rc%nt)
  real, intent(inout) :: maskarray(LDT_rc%lnc(n),LDT_rc%lnr(n))
!
! !DESCRIPTION:
!  This subroutine reads the SACHTET-v3.5.6 landcover data and returns the 
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
   integer  :: c,r,nr,nc
   integer  :: glpnc, glpnr             ! Parameter (global) total columns and rows
   integer  :: rc
   logical  :: file_exists
   integer  :: file_status
   real, allocatable  :: data2d(:,:)
   real, allocatable  :: lat(:,:)
   real, allocatable  :: lon(:,:)
   integer :: nrow, ncol, row, col
   integer :: subrow, subcol
   real    :: rlon, rlat

   real    :: vegtype(LDT_rc%lnc(n),LDT_rc%lnr(n))
   real    :: vegcnt(LDT_rc%lnc(n),LDT_rc%lnr(n),num_types)

!__________________________________________________________________

!- Assign additional land cover types, including generic water points: 
   select case ( LDT_rc%lc_type(n) )
    case ( "UMD" )
      LDT_rc%bareclass    = 12
      LDT_rc%urbanclass   = 13
      LDT_rc%waterclass   = 14
      LDT_rc%snowclass    = 0
      LDT_rc%wetlandclass = 0
      LDT_rc%glacierclass = 0
    case default ! non-supported options
      write(LDT_logunit,*) "The land classification: ",trim(LDT_rc%lc_type(n)),&
                           " does not exist for SAC-HTET 3.5.6 source."
      write(LDT_logunit,*) " -- Please select either: UMD or IGBP "
      write(LDT_logunit,*) "Program stopping ..."
      call LDT_endrun
   end select

!- Check if land cover file exists:
   inquire( file=trim(LDT_rc%vfile(n)), exist=file_exists )
   if(.not. file_exists) then
      write(LDT_logunit,*) "Landcover map: ",trim(LDT_rc%vfile(n))," does not exist "
      write(LDT_logunit,*) "Program stopping ..."
      call LDT_endrun
   endif

   inquire( file=trim(LDT_rc%mfile(n)), exist=file_exists )
   if(.not. file_exists) then
      write(LDT_logunit,*) "Landmask map: ",trim(LDT_rc%mfile(n))," does not exist "
      write(LDT_logunit,*) "Program stopping ..."
      call LDT_endrun
   endif

   write(LDT_logunit,*) "Reading SAC-HTET 3.5.6 Veg file: "//trim(LDT_rc%vfile(n))

   vegtype     = float(LDT_rc%waterclass)
   vegcnt      = 0.
   fgrd        = 0.
   maskarray   = 0.

   call LDT_transform_xmrgparam( n, LDT_rc%lnc(n), LDT_rc%lnr(n),     &
            LDT_rc%gridDesc(n,:), LDT_rc%vfile(n), LDT_rc%udef, vegtype )
            
   write(LDT_logunit,*) "Done reading: "//trim(LDT_rc%vfile(n))

!- Estimate fraction of grid (fgrid) represented by vegetation type::
!- Bring 2-D Vegtype to 3-D Vegcnt tile space:
   if ( LDT_rc%lc_gridtransform(n) == "none"     .or. &  ! -- NON-TILED SURFACES
        LDT_rc%lc_gridtransform(n) == "neighbor" .or. & 
        LDT_rc%lc_gridtransform(n) == "mode" ) then  
     do r = 1, LDT_rc%lnr(n)
        do c = 1, LDT_rc%lnc(n)
           if ( vegtype(c,r) .le. 0 ) then
              vegtype(c,r) = float(LDT_rc%waterclass)
           endif
           if( (nint(vegtype(c,r)) .ne. LDT_rc%waterclass ) .and. &
               (nint(vegtype(c,r)) .ne. LDT_rc%udef)) then
              vegcnt(c,r,NINT(vegtype(c,r))) = 1.0
           endif
        enddo
     end do
   endif   ! End NON-TILED vegetation option

   call param_index_fgrdcalc( n, LDT_rc%lc_proj, LDT_rc%lc_gridtransform(n), &
                              LDT_rc%waterclass, num_types, vegcnt, fgrd )
! .....

   write(LDT_logunit,*) "Reading SAC-HTET 3.5.6 Mask file: "//trim(LDT_rc%mfile(n))

   call LDT_transform_xmrgparam( n, LDT_rc%lnc(n), LDT_rc%lnr(n),     &
            LDT_rc%gridDesc(n,:), LDT_rc%mfile(n), 0., maskarray )

   do r = 1, LDT_rc%lnr(n)
      do c = 1, LDT_rc%lnc(n)
      !- Rearrange mask values for universal LIS-mask values (0=water,1=land):
         if( maskarray(c,r) == -1. ) maskarray(c,r) = 0.
         if( maskarray(c,r) ==  2. ) maskarray(c,r) = 1.

      !- Ensure consistency between landmask and landcover water points:
         if((fgrd(c,r,LDT_rc%nt) > LDT_rc%gridcell_water_frac(n)) .and. &
             maskarray(c,r) == 1 ) then
            maskarray(c,r) = 0.
         endif
      end do
   end do

   write(LDT_logunit,*) "Done reading: "//trim(LDT_rc%mfile(n))

end subroutine read_SACHTET356_lc
