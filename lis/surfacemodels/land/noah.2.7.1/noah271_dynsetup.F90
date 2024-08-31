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
! !ROUTINE: noah271_dynsetup
! \label{noah271_dynsetup}
!
! !REVISION HISTORY:
!  28 Jan 2002: Sujay Kumar, Initial Specification
!
! !INTERFACE:
subroutine noah271_dynsetup(n)
! !USES: 
  use LIS_coreMod, only : LIS_rc
  use LIS_snowMod, only : LIS_snow_struc
  use noah271_lsmMod, only : noah271_struc
! 
! !DESCRIPTION: 
!  This routine sets up the time-dependent variables in Noah2.7.1
!  Currently this routine is a placeholder
! 
!EOP     
  implicit none
  integer, intent(in) :: n 

  integer   :: t

  if(LIS_rc%snowsrc(n).gt.0) then 
     do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
        LIS_snow_struc(n)%sneqv(t)     = noah271_struc(n)%noah(t)%sneqv
        LIS_snow_struc(n)%snowdepth(t) = noah271_struc(n)%noah(t)%snowh
     enddo
  endif
  return
end subroutine noah271_dynsetup
