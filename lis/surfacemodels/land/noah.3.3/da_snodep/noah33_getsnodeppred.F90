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
! !ROUTINE: noah33_getsnodeppred
! \label{noah33_getsnodeppred}
!
! !REVISION HISTORY:
! 27Feb2005: Sujay Kumar; Initial Specification
! 25Jun2006: Sujay Kumar: Updated for the ESMF design
! 02 Mar 2010: Sujay Kumar; Modified for Noah 3.1
! 01 May 2014: Yuqiong Liu; modifed to include mesh8, mesh16, and 0p25 SNODEP data
!
! !INTERFACE:
subroutine noah33_getsnodeppred(n, obs_pred)

! !USES:
  use ESMF
  use LIS_coreMod, only : LIS_rc,LIS_surface
  use noah33_lsmMod
  use SNODEPobs_Mod, only: SNODEP_obs_obj

  implicit none
! !ARGUMENTS: 
  integer, intent(in)    :: n
  real                   :: obs_pred(LIS_rc%ngrid(n),LIS_rc%nensem(n))
!EOP

  integer                :: i,t,m,gid

  do i=1,LIS_rc%npatch(n,LIS_rc%lsm_index),LIS_rc%nensem(n)
     do m=1,LIS_rc%nensem(n)
        t = i+m-1
        gid = LIS_surface(n,1)%tile(t)%index
        if (SNODEP_obs_obj(n)%mesh .eq. 8) then
           obs_pred(gid,m)= noah33_struc(n)%noah(t)%snowh*39.37 !convert from meter to inch
                 !mesh16 and 0p25 deg data are in meters so no conversion needed
        elseif (SNODEP_obs_obj(n)%mesh .eq. 16 .or. SNODEP_obs_obj(n)%mesh .eq. 25) then
           obs_pred(gid,m)= noah33_struc(n)%noah(t)%snowh
        endif
     enddo
  enddo
  
end subroutine noah33_getsnodeppred

