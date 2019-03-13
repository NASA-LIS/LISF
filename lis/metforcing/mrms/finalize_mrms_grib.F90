!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.0
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LIS_misc.h"
!BOP
! !MODULE: finalize_mrms
!  \label{finalize_mrms}
!
! !REVISION HISTORY: 
! 25Oct2005; Sujay Kumar, Initial Code
! 17Jul2006; K. Arsenault, Added Stage IV
! 13 Feb 2015; Jonathan Case, modified for MRMS QPE
! 07 Sep 2017; Jessica Erlingis, modified for operational MRMS QPE
! 22 Feb 2019; Jessica Erlingis, modified to add nearest neighbor case
!
! !INTERFACE:
subroutine finalize_mrms_grib(findex)

! !USES:
  use LIS_coreMod, only : LIS_rc
  use mrms_grib_forcingMod
!
! !DESCRIPTION:
!  Routine to cleanup MRMS forcing related memory allocations.   
! 
!EOP
  implicit none
  integer, intent(IN) :: findex
  integer   :: n

! J.Case -- temporary variables
  character*50 :: dom
  real         :: res
  
  do n=1,LIS_rc%nnest
    dom=LIS_rc%lis_map_proj
    res=LIS_rc%gridDesc(n,9)
    if ( (((dom.eq."latlon").or.(dom.eq."gaussian")) .and. (res.eq.0.01)) .or. & !Nearest neighbor case
      (((dom.eq."mercator").or.(dom.eq."lambert").or.                       &
      (dom.eq."polar").or.(dom.eq."UTM")) .and. (res.eq.1.0)) ) then
      deallocate(mrms_grib_struc(n)%n113)
    endif 
    if ( (((dom.eq."latlon").or.(dom.eq."gaussian")) .and. (res.lt.0.01)) .or. &
         (((dom.eq."mercator").or.(dom.eq."lambert").or.                       &
           (dom.eq."polar").or.(dom.eq."UTM")) .and. (res.lt.1.0)) ) then
      if(trim(LIS_rc%met_interp(findex)).eq."bilinear") then 
         deallocate(mrms_grib_struc(n)%rlat1)
         deallocate(mrms_grib_struc(n)%rlon1)
         deallocate(mrms_grib_struc(n)%n111)
         deallocate(mrms_grib_struc(n)%n121)
         deallocate(mrms_grib_struc(n)%n211)
         deallocate(mrms_grib_struc(n)%n221)
         deallocate(mrms_grib_struc(n)%w111)
         deallocate(mrms_grib_struc(n)%w121)
         deallocate(mrms_grib_struc(n)%w211)
         deallocate(mrms_grib_struc(n)%w221)
      elseif(trim(LIS_rc%met_interp(findex)).eq."budget-bilinear") then 
         deallocate(mrms_grib_struc(n)%rlat2)
         deallocate(mrms_grib_struc(n)%rlon2)
         deallocate(mrms_grib_struc(n)%n112)
         deallocate(mrms_grib_struc(n)%n122)
         deallocate(mrms_grib_struc(n)%n212)
         deallocate(mrms_grib_struc(n)%n222)
         deallocate(mrms_grib_struc(n)%w112)
         deallocate(mrms_grib_struc(n)%w122)
         deallocate(mrms_grib_struc(n)%w212)
         deallocate(mrms_grib_struc(n)%w222)
      endif
    else ! upscaling deallocate
      deallocate(mrms_grib_struc(n)%n11)
    endif
 enddo
 deallocate(mrms_grib_struc)

end subroutine finalize_mrms_grib
