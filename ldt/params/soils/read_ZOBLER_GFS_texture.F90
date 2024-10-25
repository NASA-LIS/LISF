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
! !ROUTINE: read_ZOBLER_GFS_texture
! \label{read_ZOBLER_GFS_texture}
!
! !REVISION HISTORY:
!  May 2014: Grey Nearing; Initial Specification
!
! !INTERFACE:
subroutine read_ZOBLER_GFS_texture(n,num_bins,array,texture_layers)

! !USES:
  use LDT_coreMod,    only : LDT_rc
  use LDT_logMod,     only : LDT_logunit, LDT_getNextUnitNumber, &
       LDT_releaseUnitNumber, LDT_endrun
  use LDT_fileIOMod,  only : readLISdata, LDT_checkDomainExtents 
  use LDT_paramTileInputMod, only: param_index_fgrdcalc
  use LDT_gridmappingMod

  implicit none
! !ARGUMENTS: 
  integer,intent(in)    :: n
  integer,intent(in)    :: num_bins   ! Number of soil types
  real,   intent(inout) :: array(LDT_rc%lnc(n),LDT_rc%lnr(n),num_bins)
  real,   intent(inout) :: texture_layers(LDT_rc%lnc(n),LDT_rc%lnr(n),1)
!
! !DESCRIPTION:
!  This subroutine retrieves GFS soil texture data and reprojects
!  it to the latlon projection. 
!
!  The arguments are:
!  \begin{description}
!  \item[n]
!   index of the nest
!  \item[array]
!   output field with the retrieved soil texture
!  \end{description}
!EOP
  integer :: c, r, t 
  integer :: water_class,isum
  integer :: ftn
  !real    :: fgrd(LDT_rc%lnc(n),LDT_rc%lnr(n),num_bins)
  real, allocatable :: fgrd(:,:,:)
  logical :: file_exists
  integer :: glpnc, glpnr             ! Parameter (global) total columns and rows
  integer :: subpnc, subpnr           ! Parameter subsetted columns and rows
  real    :: subparam_gridDesc(20)    ! Input parameter grid desc array
  integer, allocatable  :: lat_line(:,:), lon_line(:,:)
  real,    allocatable  :: dummy_lat(:,:), dummy_lon(:,:), read_soil(:,:)
! ___________________________________________________________________

  allocate(fgrd(LDT_rc%lnc(n),LDT_rc%lnr(n),num_bins))
  fgrd = 0.
  array = 0.
  texture_layers = 0.

  inquire(file=trim(LDT_rc%txtfile(n)), exist=file_exists)
  if(.not.file_exists) then 
     write(LDT_logunit,*) "Texture map ",trim(LDT_rc%txtfile(n))," not found."
     write(LDT_logunit,*) "Program stopping ..."
     call LDT_endrun
  endif

  write(unit=LDT_logunit,fmt=*) "[INFO] Reading ZOBLER GFS texture file: ",&
                                trim(LDT_rc%txtfile(n))

! - find subparam grid
   subparam_gridDesc = 0.
   call LDT_RunDomainPts(n,LDT_rc%soiltext_proj,LDT_rc%soiltext_gridDesc(n,:), &
             glpnc,glpnr,subpnc,subpnr,subparam_gridDesc,lat_line,lon_line)

!- Open file:
   ftn = LDT_getNextUnitNumber()
   open(ftn,file=LDT_rc%txtfile(n),form='unformatted',status='old',&
        access='sequential',recl=4)

! - read file
   allocate(dummy_lat(glpnc,glpnr))
   allocate(dummy_lon(glpnc,glpnr))
   allocate(read_soil(glpnc,glpnr))
   read(ftn) dummy_lat
   read(ftn) dummy_lon
   read(ftn) read_soil
   close(ftn)
   deallocate(dummy_lat)
   deallocate(dummy_lon)

!- Transform parameter grid to LIS run domain:
   select case (LDT_rc%soiltext_gridtransform(n))
     case("none","tile")
       do r = 1,subpnr
         do c = 1,subpnc
           do t = 1,num_bins
             if (read_soil(lon_line(c,r),lat_line(c,r)).eq.t) then
               fgrd(c,r,t) = 1.
             endif
           enddo ! tilespace
         enddo ! columns
       enddo ! rows
       deallocate(read_soil)
    case default
       write(LDT_logunit,*) " This spatial transform, ",trim(LDT_rc%soiltext_gridtransform(n)),&
           ", for GFS soil texture is not available at this time ... "
      write(LDT_logunit,*) "Program stopping ..."
      call LDT_endrun
   end select
   close(ftn)

! - add water class
!   do r = 1,subpnr
!     do c = 1,subpnc
!       isum = sum(fgrd(c,r,:))
!       fgrd(c,r,water_class) = 1.-isum
!     enddo
!   enddo

!- Estimate fraction of grid (fgrid) represented by soil type::
   array = fgrd
   deallocate(fgrd)

! ---
   call LDT_releaseUnitNumber(ftn)
   write(unit=LDT_logunit,fmt=*) "[INFO] Done reading ZOBLER GFS soil texture file."

 end subroutine read_ZOBLER_GFS_texture

