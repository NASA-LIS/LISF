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
! !ROUTINE: noahmp401_getsnowpred
! \label{noahmp401_getsnowpred}
!
! !REVISION HISTORY:
! 27Feb2005: Sujay Kumar; Initial Specification
! 25Jun2006: Sujay Kumar: Updated for the ESMF design
! 02 Mar 2010: Sujay Kumar; Modified for Noah 3.1
! 01 May 2014: Yuqiong Liu; modifed to include mesh8, mesh16, and 0p25 SNODEP data
! 24 May 2017: Yeosang Yoon: updated the file to work with the DA observation
!              space updates.
! 03 Oct 2018: Yeosang Yoon; Modified for NoahMP 3.6
! 14 Dec 2018: Yeosang Yoon; Modified for NoahMP 4.0.1 and SNODEP
! 15 May 2019: Yeosang Yoon; Modified for NoahMP 4.0.1 and LDTSI
! 13 Dec 2019: Eric Kemp; Replaced LDTSI with SNOW
!
! !INTERFACE:
subroutine noahmp401_getsnowpred(n, k, obs_pred)

! !USES:
  use ESMF
  use LIS_coreMod, only : LIS_rc,LIS_surface
  use noahmp401_lsmMod
  use LIS_DAobservationsMod

  implicit none
! !ARGUMENTS: 
  integer, intent(in)    :: n
  integer, intent(in)    :: k
  real                   :: obs_pred(LIS_rc%ngrid(n),LIS_rc%nensem(n))
  real                   :: snwd(LIS_rc%npatch(n,LIS_rc%lsm_index))
!EOP

! !DESCRIPTION:
!  This routine computes the obspred ('Hx') term for SNOW DA assimilation
!  instances.

  integer                :: t

  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     snwd(t) = noahmp401_struc(n)%noahmp401(t)%snowh ! Keep in meters
  enddo

  call LIS_convertPatchSpaceToObsEnsSpace(n,k,&
       LIS_rc%lsm_index, &
       snwd,&
       obs_pred)
  
end subroutine noahmp401_getsnowpred

