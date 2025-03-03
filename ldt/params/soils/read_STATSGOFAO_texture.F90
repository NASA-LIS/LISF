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
! !ROUTINE: read_STATSGOFAO_texture
! \label{read_STATSGOFAO_texture}
!
! !REVISION HISTORY:
!  03 Sept 2004: Sujay Kumar; Initial Specification
!
! !INTERFACE:
subroutine read_STATSGOFAO_texture( n, num_bins, array, texture_layers )

! !USES:
  use LDT_coreMod,    only : LDT_rc
  use LDT_logMod,     only : LDT_logunit, LDT_getNextUnitNumber, &
       LDT_releaseUnitNumber, LDT_endrun
  use LDT_fileIOMod,  only : readLISdata, LDT_checkDomainExtents 
  use LDT_paramTileInputMod, only: param_index_fgrdcalc

  implicit none
! !ARGUMENTS: 
  integer,intent(in)    :: n
  integer,intent(in)    :: num_bins   ! Number of soil types
  real,   intent(inout) :: array(LDT_rc%lnc(n),LDT_rc%lnr(n),num_bins)
  real,   intent(inout) :: texture_layers(LDT_rc%lnc(n),LDT_rc%lnr(n),1)
!
! !DESCRIPTION:
!  This subroutine retrieves STATSGO soil texture data and reprojects
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
  integer :: water_class
  integer :: ftn
  real    :: temp(LDT_rc%lnc(n),LDT_rc%lnr(n),1)
  real    :: fgrd(LDT_rc%lnc(n),LDT_rc%lnr(n),num_bins)
  logical :: file_exists
! ___________________________________________________________________
  
  water_class = 14    ! Water class for NCAR STATSGO/FAO soil texture class

  temp = LDT_rc%udef
  fgrd = 0.
  array = 0.
  texture_layers = 0.

  inquire(file=trim(LDT_rc%txtfile(n)), exist=file_exists)
  if(.not.file_exists) then 
     write(LDT_logunit,*) "Texture map ",trim(LDT_rc%txtfile(n))," not found."
     write(LDT_logunit,*) "Program stopping ..."
     call LDT_endrun
  endif
  write(unit=LDT_logunit,fmt=*) "[INFO] Reading STATSGO texture file: ",&
                                trim(LDT_rc%txtfile(n))
!- Open file:
   ftn = LDT_getNextUnitNumber()
   open( ftn,file=LDT_rc%txtfile(n),form='unformatted',status='old',&
        access='direct',recl=4 )

!- Transform parameter grid to LIS run domain:
   select case ( LDT_rc%soiltext_gridtransform(n) )

 !- Selecting only dominant soil types for 3D array: 
    case( "none", "mode", "neighbor" )
      call readLISdata(n, ftn, LDT_rc%soiltext_proj, &
           LDT_rc%soiltext_gridtransform(n), &
           LDT_rc%soiltext_gridDesc(n,:), 1, temp )  ! 1 indicates 2D layer
      
   !- Bring 2-D Array to 3-D Soil tile space:
      do r = 1, LDT_rc%lnr(n)
         do c = 1, LDT_rc%lnc(n)
         !- Assign soil texture types of less than 0 to an actual texture value:
            if( nint(temp(c,r,1)) .le. 0 ) then
              temp(c,r,1) = 12   ! STATSGO -- Clay 
            endif
            if( (nint(temp(c,r,1)) .ne. water_class  ) .and. &
                (nint(temp(c,r,1)) .ne. LDT_rc%udef) ) then
               array(c,r,NINT(temp(c,r,1))) = 1.0
            endif
         enddo
      enddo

 !- Determining soil texture fraction of each gridcell for 3D array: 
    case( "tile" )
      call readLISdata(n, ftn, LDT_rc%soiltext_proj, LDT_rc%soiltext_gridtransform(n), &
                       LDT_rc%soiltext_gridDesc(n,:), num_bins, array )

    case default
      write(LDT_logunit,*) " This spatial transform, ",trim(LDT_rc%soiltext_gridtransform(n)),&
                           " , for STATSGO+FAO soil texture is not available at this time ... "
      write(LDT_logunit,*) "Program stopping ..."
      call LDT_endrun

   end select
   close(ftn)
  
!- Estimate fraction of grid (fgrid) represented by soil type::
   call param_index_fgrdcalc( n, LDT_rc%soiltext_proj, LDT_rc%soiltext_gridtransform(n), &
                              water_class, num_bins, array, fgrd )
   array = fgrd

! ---
   call LDT_releaseUnitNumber(ftn)
   write(unit=LDT_logunit,fmt=*) "[INFO] Done reading STATSGO+FAO soil texture file."

 end subroutine read_STATSGOFAO_texture

