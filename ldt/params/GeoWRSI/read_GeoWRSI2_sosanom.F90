!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
!
! !ROUTINE: read_GeoWRSI2_sosanom
!  \label{read_GeoWRSI2_sosanom}

! !REVISION HISTORY:
!  25 Oct 2013: K. Arsenault; Initial Specification
!
! !INTERFACE:
subroutine read_GeoWRSI2_sosanom(n,array)

! !USES:
  use ESMF
  use LDT_coreMod, only : LDT_rc, LDT_domain
  use LDT_logMod,  only : LDT_logunit, LDT_getNextUnitNumber, &
          LDT_releaseUnitNumber, LDT_endrun
  use map_utils
  use geowrsi2_arraymgmt_module, only: nullify_ptr, alloc_arr, dealloc_arr
  use fbil_module
  use GeoWRSI_parmsMod
!EOP      
  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: n
  real, intent(inout) :: array(LDT_rc%lnc(n),LDT_rc%lnr(n))
!
! !DESCRIPTION:
!  This subroutine retrieves the WRSI model SOS anomaly
!   (SOSanom) parameter from a *BIL formatted set of files (includes
!   *hdr file). 
!
!  The arguments are:
!  \begin{description}
!   \item[n]
!    index of the nest
!   \item[array]
!    output field for SOSanom data
!  \end{description}
!
!EOP      
  integer  :: ftn, rc
  integer  :: nc, nr
  logical  :: file_exists
  logical  :: fileread_ok
  logical  :: real4ptr2dIsNull

  real*4,    pointer :: real4ptr2d(:,:)
  type(charN)        :: bil_filename

! _________________________________________________________________________

   inquire( file=trim(GeoWRSI_struc(n)%sosanom_file)//".bil", exist=file_exists )
   if(.not. file_exists) then
      write(LDT_logunit,*) "The WRSI SOSanom file, ",trim(GeoWRSI_struc(n)%sosanom_file),&
                           " does not exist."
      write(LDT_logunit,*) "Program stopping ..."
      call LDT_endrun
   endif

!- Setup BIL Geographic Coordination information (from header files):
!  Allocate (internal I/O lib code) keeper of the model grid information
   allocate(gCoords)

 ! Set the grid space
   call geowrsi2_set_gcoords(n)

   call calcXyBounds( &
                      addOffset_arg=(0+1),      &
                      ifMidPointRoundUp=.true., &
                      minLon=gCoords%minLon,    &
                      ulxmap=gCoords%minLon,    &
                      xdim=gCoords%pixLon,      &
                      maxLon=gCoords%maxLon,    &
                      ulymap=gCoords%maxLat,    &
                      minLat=gCoords%minLat,    &
                      ydim=gCoords%pixLat,      &
                      maxLat=gCoords%maxLat,    &
                      minX=gCoords%minX,        &
                      maxX=gCoords%maxX,        &
                      minY=gCoords%minY,        &
                      maxY=gCoords%maxY         )


 ! Allocate temporary memory for staging the data in:
   call nullify_ptr(real4ptr2dIsNull, real4ptr2d=real4ptr2d)

   real4ptr2dIsNull = alloc_arr( real4ptr2dIsNull, &
                                 dim1Sz=((gCoords%maxX-gCoords%minX)+1), &
                                 dim2Sz=((gCoords%maxY-gCoords%minY)+1), &
                                 real4ptr2d=real4ptr2d )

!- Assign filename pointer:
   allocate(bil_filename%str)
   bil_filename%str = GeoWRSI_struc(n)%sosanom_file
   fileread_ok = populate_array_from_bilfile_2d(    &
                          mapFilename=bil_filename, &
                          real4ptr2d=real4ptr2d )

   if( fileread_ok ) then
      write(LDT_logunit,*) "[INFO] WRSI-Parameter: SOS-anomaly read-in successfully."
   else
      write(LDT_logunit,*) "[ERR] SOS-anomaly NOT read-in."
      write(LDT_logunit,*) "Stopping ..."; call LDT_endrun
      real4ptr2d = LDT_rc%udef
   endif

   array = LDT_rc%udef
   do nr = 1, LDT_rc%lnr(n)
      do nc = 1, LDT_rc%lnc(n)
         array(nc,nr) = real4ptr2d(nc,nr)
      end do
   end do

! --- 
   call dealloc_arr(real4ptr2dIsNull, real4ptr2d=real4ptr2d)
   call nullify_ptr(real4ptr2dIsNull, real4ptr2d=real4ptr2d)
   deallocate(bil_filename%str)
   deallocate(gCoords)

end subroutine read_GeoWRSI2_sosanom
