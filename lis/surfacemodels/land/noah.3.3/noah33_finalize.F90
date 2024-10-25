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
! !ROUTINE: noah33_finalize
! \label{noah33_finalize}
!
! !REVISION HISTORY:
!  28 Apr 2002: Sujay Kumar; Initial code
!   8 May 2009: Sujay Kumar; additions for Noah3.1
!  27 Oct 2010: David Mocko, changes for Noah3.1 in LIS6.1
!   7 Nov 2010: David Mocko, changes for Noah3.2 in LIS6.1
!   9 Sep 2011: David Mocko, changes for Noah3.3 in LIS6.1
!
! !INTERFACE:
subroutine noah33_finalize()
! !USES:
  use LIS_coreMod,  only : LIS_rc
  use noah33_lsmMod
!
! !DESCRIPTION:
!  
!  This routine cleans up the allocated memory structures in Noah3.3
!  
!EOP
  implicit none

  integer :: t,n

  do n=1,LIS_rc%nnest
     do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
        deallocate(noah33_struc(n)%noah(t)%stc    )
        deallocate(noah33_struc(n)%noah(t)%smc    )
        deallocate(noah33_struc(n)%noah(t)%sh2o   )  
        deallocate(noah33_struc(n)%noah(t)%relsmc )  
     enddo
     deallocate(noah33_struc(n)%noah)
     deallocate(noah33_struc(n)%lyrthk)
     deallocate(noah33_struc(n)%inittemp)
     deallocate(noah33_struc(n)%initsm)
     deallocate(noah33_struc(n)%initsmliq)
     deallocate(noah33_struc(n)%shdtbl   )
     deallocate(noah33_struc(n)%nrotbl   )
     deallocate(noah33_struc(n)%rstbl    )
     deallocate(noah33_struc(n)%rgltbl   )
     deallocate(noah33_struc(n)%hstbl    )
     deallocate(noah33_struc(n)%snuptbl  )
     deallocate(noah33_struc(n)%maxalb   )
     deallocate(noah33_struc(n)%emissmin )
     deallocate(noah33_struc(n)%emissmax )
     deallocate(noah33_struc(n)%laimin   )
     deallocate(noah33_struc(n)%laimax   )
     deallocate(noah33_struc(n)%albmin   )
     deallocate(noah33_struc(n)%albmax   )
     deallocate(noah33_struc(n)%z0min    )
     deallocate(noah33_struc(n)%z0max    ) 
  enddo
  deallocate(noah33_struc)  
end subroutine noah33_finalize
