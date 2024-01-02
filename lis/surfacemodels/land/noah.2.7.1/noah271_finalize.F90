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
! !ROUTINE: noah271_finalize
! \label{noah271_finalize}
!
! !REVISION HISTORY:
!  28 Apr 2002: Sujay Kumar; Initial code
!  27 Oct 2010: David Mocko, changes for Noah2.7.1 in LIS6.1
! 
! !INTERFACE:
subroutine noah271_finalize()
! !USES:
  use LIS_coreMod,  only : LIS_rc
  use noah271_lsmMod
!
! !DESCRIPTION:
!  
!  This routine cleans up the allocated memory structures in Noah2.7.1
!  
!EOP
  implicit none
  integer :: t,n
  
  do n=1,LIS_rc%nnest
     do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
        deallocate(noah271_struc(n)%noah(t)%stc)
        deallocate(noah271_struc(n)%noah(t)%smc)
        deallocate(noah271_struc(n)%noah(t)%sh2o)
        deallocate(noah271_struc(n)%noah(t)%relsmc )  
     enddo
     deallocate(noah271_struc(n)%lyrthk)
     deallocate(noah271_struc(n)%inittemp)
     deallocate(noah271_struc(n)%initsm)
     deallocate(noah271_struc(n)%initsmliq)
     deallocate(noah271_struc(n)%noah)
  enddo
  deallocate(noah271_struc)

end subroutine noah271_finalize
