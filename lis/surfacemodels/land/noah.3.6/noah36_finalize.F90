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
! !ROUTINE: noah36_finalize
! \label{noah36_finalize}
!
! !REVISION HISTORY:
!  28 Apr 2002: Sujay Kumar; Initial code
!   8 May 2009: Sujay Kumar; additions for Noah3.1
!  27 Oct 2010: David Mocko, changes for Noah3.1 in LIS6.1
!   7 Nov 2010: David Mocko, changes for Noah3.2 in LIS6.1
!   9 Sep 2011: David Mocko, changes for Noah3.3 in LIS6.1
!  30 Oct 2014: David Mocko, added Noah-3.6 into LIS-7
!
! !INTERFACE:
subroutine noah36_finalize()
! !USES:
  use LIS_coreMod,  only : LIS_rc
  use noah36_lsmMod
!
! !DESCRIPTION:
!  
!  This routine cleans up the allocated memory structures in Noah-3.6
!  
!EOP
  implicit none

  integer :: t,n

  do n=1,LIS_rc%nnest
     do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
        deallocate(noah36_struc(n)%noah(t)%stc    )
        deallocate(noah36_struc(n)%noah(t)%smc    )
        deallocate(noah36_struc(n)%noah(t)%sh2o   )  
        deallocate(noah36_struc(n)%noah(t)%relsmc )  
     enddo
     deallocate(noah36_struc(n)%noah)
     deallocate(noah36_struc(n)%lyrthk)
     deallocate(noah36_struc(n)%inittemp)
     deallocate(noah36_struc(n)%initsm)
     deallocate(noah36_struc(n)%initsmliq)
     deallocate(noah36_struc(n)%shdtbl   )
     deallocate(noah36_struc(n)%nrotbl   )
     deallocate(noah36_struc(n)%rstbl    )
     deallocate(noah36_struc(n)%rgltbl   )
     deallocate(noah36_struc(n)%hstbl    )
     deallocate(noah36_struc(n)%snuptbl  )
     deallocate(noah36_struc(n)%maxalb   )
     deallocate(noah36_struc(n)%emissmin )
     deallocate(noah36_struc(n)%emissmax )
     deallocate(noah36_struc(n)%laimin   )
     deallocate(noah36_struc(n)%laimax   )
     deallocate(noah36_struc(n)%albmin   )
     deallocate(noah36_struc(n)%albmax   )
     deallocate(noah36_struc(n)%z0min    )
     deallocate(noah36_struc(n)%z0max    ) 
  enddo
  deallocate(noah36_struc)  
end subroutine noah36_finalize
