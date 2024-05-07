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
! !MODULE: finalize_gfs
! \label{gfsforcing_finalie}
! 
! !REVISION HISTORY: 
!  16 Mar 2008: Sujay Kumar; Initial specification
! 
! !INTERFACE:
subroutine finalize_gfs(findex)
! !USES:
  use LIS_coreMod,  only : LIS_rc
  use gfs_forcingMod
!
! !DESCRIPTION:
!  Routine to cleanup allocated structures for GFS forcing. 
!
!EOP  
  implicit none
  integer :: findex
  integer :: n 
  
  do n=1,LIS_rc%nnest  
     if(trim(LIS_rc%met_interp(findex)).eq."bilinear") then 

        deallocate(gfs_struc(n)%n111)
        deallocate(gfs_struc(n)%n121)
        deallocate(gfs_struc(n)%n211)
        deallocate(gfs_struc(n)%n221)
        deallocate(gfs_struc(n)%w111)
        deallocate(gfs_struc(n)%w121)
        deallocate(gfs_struc(n)%w211)
        deallocate(gfs_struc(n)%w221)
     elseif(trim(LIS_rc%met_interp(findex)).eq."budget-bilinear") then 

        deallocate(gfs_struc(n)%n111)
        deallocate(gfs_struc(n)%n121)
        deallocate(gfs_struc(n)%n211)
        deallocate(gfs_struc(n)%n221)
        deallocate(gfs_struc(n)%w111)
        deallocate(gfs_struc(n)%w121)
        deallocate(gfs_struc(n)%w211)
        deallocate(gfs_struc(n)%w221)

        deallocate(gfs_struc(n)%n112)
        deallocate(gfs_struc(n)%n122)
        deallocate(gfs_struc(n)%n212)
        deallocate(gfs_struc(n)%n222)
        deallocate(gfs_struc(n)%w112)
        deallocate(gfs_struc(n)%w122)
        deallocate(gfs_struc(n)%w212)
        deallocate(gfs_struc(n)%w222)
     endif
  enddo
  deallocate(gfs_struc)
end subroutine finalize_gfs
