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
! !ROUTINE: finalize_agrradps
! \label{finalize_agrradps}
!
! !REVISION HISTORY: 
! 25Oct2005; Sujay Kumar, Initial Code
! 
! !INTERFACE:
subroutine finalize_agrradps
! !USES:
  use LIS_coreMod,         only : LIS_rc
  use agrradps_forcingMod, only : agrradps_struc
!
! !DESCRIPTION:
!  Routine to cleanup AGRRADPS forcing related memory allocations.   
! 
!EOP
  implicit none
  
  integer   :: n
  
  do n=1,LIS_rc%nnest
     deallocate(agrradps_struc(n)%smask)
     if(agrradps_struc(n)%gridspan.eq.1) then
       deallocate(agrradps_struc(n)%rlat1_nh)
       deallocate(agrradps_struc(n)%rlon1_nh)
       deallocate(agrradps_struc(n)%n111_nh)
       deallocate(agrradps_struc(n)%n121_nh)
       deallocate(agrradps_struc(n)%n211_nh)
       deallocate(agrradps_struc(n)%n221_nh)
       deallocate(agrradps_struc(n)%w111_nh)
       deallocate(agrradps_struc(n)%w121_nh)
       deallocate(agrradps_struc(n)%w211_nh)
       deallocate(agrradps_struc(n)%w221_nh)
     elseif(agrradps_struc(n)%gridspan.eq.2) then
       deallocate(agrradps_struc(n)%rlat1_sh)
       deallocate(agrradps_struc(n)%rlon1_sh)
       deallocate(agrradps_struc(n)%n111_sh)
       deallocate(agrradps_struc(n)%n121_sh)
       deallocate(agrradps_struc(n)%n211_sh)
       deallocate(agrradps_struc(n)%n221_sh)
       deallocate(agrradps_struc(n)%w111_sh)
       deallocate(agrradps_struc(n)%w121_sh)
       deallocate(agrradps_struc(n)%w211_sh)
       deallocate(agrradps_struc(n)%w221_sh)
     endif
 enddo
 deallocate(agrradps_struc)

end subroutine finalize_agrradps
