!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
! !ROUTINE: noah39_getldtsipred
! \label{noah39_getldtsipred}
!
! !REVISION HISTORY:
! 27Feb2005: Sujay Kumar; Initial Specification
! 25Jun2006: Sujay Kumar: Updated for the ESMF design
! 02 Mar 2010: Sujay Kumar; Modified for Noah 3.1
! 01 May 2014: Yuqiong Liu; modifed to include mesh8, mesh16, and 0p25 SNODEP data
! 24 May 2017: Yeosang Yoon: updated the file to work with the DA observation
!              space updates. 
! 09 Apr 2019: Eric Kemp: Updated for Noah 3.9 and LDT-SI
!
! !INTERFACE:
subroutine noah39_getldtsipred(n, k, obs_pred)

! !USES:
  use ESMF
  use LIS_coreMod, only : LIS_rc
  use noah39_lsmMod
  use LIS_DAobservationsMod

  implicit none
! !ARGUMENTS: 
  integer, intent(in)    :: n
  integer, intent(in)    :: k
  real                   :: obs_pred(LIS_rc%ngrid(n),LIS_rc%nensem(n))
!
! !DESCRIPTION: 
!  This routine computes the obspred ('Hx') term for LDT-SI DA assimilation
!  instances. 
! 
!EOP

  real                   :: snwd(LIS_rc%npatch(n,LIS_rc%lsm_index))

  integer                :: t

 
  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     snwd(t) = noah39_struc(n)%noah(t)%snowh ! Keep in meters
  enddo

  call LIS_convertPatchSpaceToObsEnsSpace(n,k,&
       LIS_rc%lsm_index, &
       snwd,&
       obs_pred)

end subroutine noah39_getldtsipred

