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
! !ROUTINE: noah32_finalize
! \label{noah32_finalize}
!
! !REVISION HISTORY:
!  28 Apr 2002: Sujay Kumar; Initial code
!   8 May 2009: Sujay Kumar; additions for Noah3.1
!  27 Oct 2010: David Mocko, changes for Noah3.1 in LIS6.1
!   7 Nov 2010: David Mocko, changes for Noah3.2 in LIS6.1
!
! !INTERFACE:
subroutine noah32_finalize()
! !USES:
  use LIS_coreMod,  only : LIS_rc
  use noah32_lsmMod
!
! !DESCRIPTION:
!  
!  This routine cleans up the allocated memory structures in Noah3.2
!  
!EOP
  implicit none

  integer :: t,n

  do n=1,LIS_rc%nnest
     do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
        deallocate(noah32_struc(n)%noah(t)%stc)
        deallocate(noah32_struc(n)%noah(t)%smc)
        deallocate(noah32_struc(n)%noah(t)%sh2o)
        deallocate(noah32_struc(n)%noah(t)%relsmc )  
     enddo
     deallocate(noah32_struc(n)%noah)
     deallocate(noah32_struc(n)%lyrthk)
     deallocate(noah32_struc(n)%inittemp)
     deallocate(noah32_struc(n)%initsm)
     deallocate(noah32_struc(n)%initsmliq)
     deallocate(noah32_struc(n)%emissmin)
     deallocate(noah32_struc(n)%emissmax)
     deallocate(noah32_struc(n)%laimin)
     deallocate(noah32_struc(n)%laimax)
     deallocate(noah32_struc(n)%albmin)
     deallocate(noah32_struc(n)%albmax)
     deallocate(noah32_struc(n)%z0min)
     deallocate(noah32_struc(n)%z0max)
     deallocate(noah32_struc(n)%shdtbl)
     deallocate(noah32_struc(n)%nrotbl)
     deallocate(noah32_struc(n)%rstbl)
     deallocate(noah32_struc(n)%rgltbl)
     deallocate(noah32_struc(n)%hstbl)
     deallocate(noah32_struc(n)%snuptbl)
     deallocate(noah32_struc(n)%maxalb)
  enddo
  deallocate(noah32_struc)
end subroutine noah32_finalize
