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
! !ROUTINE: noah36_getsnodeppred
! \label{noah36_getsnodeppred}
!
! !REVISION HISTORY:
! 27Feb2005: Sujay Kumar; Initial Specification
! 25Jun2006: Sujay Kumar: Updated for the ESMF design
! 02 Mar 2010: Sujay Kumar; Modified for Noah 3.1
! 01 May 2014: Yuqiong Liu; modifed to include mesh8, mesh16, and 0p25 SNODEP data
! 24 May 2017: Yeosang Yoon: updated the file to work with the DA observation
!              space updates. 
!
! !INTERFACE:
subroutine noah36_getsnodeppred(n, k, obs_pred)

! !USES:
  use ESMF
  use LIS_coreMod, only : LIS_rc,LIS_surface
  use noah36_lsmMod
  use SNODEPobs_Mod, only: SNODEP_obs_obj
  use LIS_DAobservationsMod

  implicit none
! !ARGUMENTS: 
  integer, intent(in)    :: n
  integer, intent(in)    :: k
  real                   :: obs_pred(LIS_rc%ngrid(n),LIS_rc%nensem(n))
!
! !DESCRIPTION: 
!  This routine computes the obspred ('Hx') term for SNODEP DA assimilation
!  instances. The 8th mesh data is in inches whereas the 16th mesh and 
!  25km data are provided in meters. Depending on the input SNODEP 
!  data, this routine returns the obspred estimates in inches or meters. 
! 
!EOP

  real                   :: snwd(LIS_rc%npatch(n,LIS_rc%lsm_index))

  integer                :: i,t,m,gid

 
  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     if (SNODEP_obs_obj(n)%mesh .eq. 8) then
        snwd(t) = noah36_struc(n)%noah(t)%snowh*39.37 !convert from meter to inch
     elseif (SNODEP_obs_obj(n)%mesh .eq. 16 .or. SNODEP_obs_obj(n)%mesh .eq. 25) then
        snwd(t) = noah36_struc(n)%noah(t)%snowh
     endif
  enddo

  call LIS_convertPatchSpaceToObsEnsSpace(n,k,&
       LIS_rc%lsm_index, &
       snwd,&
       obs_pred)

end subroutine noah36_getsnodeppred

