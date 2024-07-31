!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
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
! 06 Jun 2019; Jessica Erlingis, modified to dellocate arrays if allocated 
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

    ! Deallocate each array if it has been previously allocated

    if ( allocated(mrms_grib_struc(n)%n113) ) then
       deallocate(mrms_grib_struc(n)%n113)
    endif
  
    if ( allocated(mrms_grib_struc(n)%rlat1) ) then
       deallocate(mrms_grib_struc(n)%rlat1)
    endif 

    if ( allocated(mrms_grib_struc(n)%rlon1) ) then
       deallocate(mrms_grib_struc(n)%rlon1)
    endif

    if ( allocated(mrms_grib_struc(n)%n111) ) then
       deallocate(mrms_grib_struc(n)%n111)
    endif

    if ( allocated(mrms_grib_struc(n)%n121) ) then
       deallocate(mrms_grib_struc(n)%n121)
    endif

    if ( allocated(mrms_grib_struc(n)%n211) ) then
       deallocate(mrms_grib_struc(n)%n211)
    endif

    if ( allocated(mrms_grib_struc(n)%n221) ) then
       deallocate(mrms_grib_struc(n)%n221)
    endif

    if ( allocated(mrms_grib_struc(n)%w111) ) then
       deallocate(mrms_grib_struc(n)%w111)
    endif

    if ( allocated(mrms_grib_struc(n)%w121) ) then
       deallocate(mrms_grib_struc(n)%w121)
    endif
    
    if ( allocated(mrms_grib_struc(n)%w211) ) then
       deallocate(mrms_grib_struc(n)%w211)
    endif

    if ( allocated(mrms_grib_struc(n)%w221) ) then
       deallocate(mrms_grib_struc(n)%w221) 
    endif

    if ( allocated(mrms_grib_struc(n)%rlat2) ) then
       deallocate(mrms_grib_struc(n)%rlat2)
    endif

    if ( allocated(mrms_grib_struc(n)%rlon2) ) then
       deallocate(mrms_grib_struc(n)%rlon2)
    endif

    if ( allocated(mrms_grib_struc(n)%n112) ) then
       deallocate(mrms_grib_struc(n)%n112)
    endif

    if ( allocated(mrms_grib_struc(n)%n122) ) then
       deallocate(mrms_grib_struc(n)%n122)
    endif

    if ( allocated(mrms_grib_struc(n)%n212) ) then
       deallocate(mrms_grib_struc(n)%n212)
    endif

    if ( allocated(mrms_grib_struc(n)%n222) ) then
       deallocate(mrms_grib_struc(n)%n222)
    endif

    if ( allocated(mrms_grib_struc(n)%w112) ) then
       deallocate(mrms_grib_struc(n)%w112)
    endif

    if ( allocated(mrms_grib_struc(n)%w122) ) then
       deallocate(mrms_grib_struc(n)%w122)
    endif

    if ( allocated(mrms_grib_struc(n)%w212) ) then 
       deallocate(mrms_grib_struc(n)%w212)
    endif
  
    if ( allocated(mrms_grib_struc(n)%w222) ) then
       deallocate(mrms_grib_struc(n)%w222)
    endif

    if ( allocated(mrms_grib_struc(n)%n11) ) then
       deallocate(mrms_grib_struc(n)%n11)
    endif

 enddo
 deallocate(mrms_grib_struc)

end subroutine finalize_mrms_grib
