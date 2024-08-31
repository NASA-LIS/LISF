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
! !ROUTINE: jules50_getusafsipred
! \label{jules50_getusafsipred}
!
! !REVISION HISTORY:
! 27Feb2005: Sujay Kumar; Initial Specification
! 25Jun2006: Sujay Kumar: Updated for the ESMF design
! 02 Mar 2010: Sujay Kumar; Modified for Noah 3.1
! 01 May 2014: Yuqiong Liu; modifed to include mesh8, mesh16, and 0p25 SNODEP data
! 24 May 2017: Yeosang Yoon: updated the file to work with the DA observation
!              space updates. 
! 05 Nov 2018: Yeosang Yoon; Modified for Jules 5.0 and SNODEP data
! 08 Jul 2019: Yeosang Yoon; Modified for Jules.5.0 and LDT-SI data
! 13 Dec 2019: Eric Kemp; Replaced LDTSI with USAFSI
!
! !INTERFACE:
subroutine jules50_getusafsipred(n, k, obs_pred)

! !USES:
  use ESMF
  use LIS_coreMod, only : LIS_rc,LIS_surface
  use jules50_lsmMod
  use LIS_DAobservationsMod

  implicit none
! !ARGUMENTS: 
  integer, intent(in)    :: n
  integer, intent(in)    :: k
  real                   :: obs_pred(LIS_rc%ngrid(n),LIS_rc%nensem(n))
!
! !DESCRIPTION: 
! This routine computes the obspred ('Hx') term for USAF DA assimilation
! instances.
! 
!EOP

  real                   :: snwd(LIS_rc%npatch(n,LIS_rc%lsm_index))
  integer                :: t, pft

  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     pft = jules50_struc(n)%jules50(t)%pft

     snwd(t) = jules50_struc(n)%jules50(t)%snowdepth(pft) ! unit: m
  enddo

  call LIS_convertPatchSpaceToObsEnsSpace(n,k,&
       LIS_rc%lsm_index, &
       snwd,&
       obs_pred)

end subroutine jules50_getusafsipred

