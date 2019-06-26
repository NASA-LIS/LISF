!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
! !ROUTINE: noah33_getsnowpred_PMWsnow
! \label{noah33_getsnowpred_PMWsnow}
!
! !REVISION HISTORY:
! 27Feb2005: Sujay Kumar; Initial Specification
! 25Jun2006: Sujay Kumar: Updated for the ESMF design
!  02 Mar 2010: Sujay Kumar; Modified for Noah 3.1
! Mar 14 2014, Yuqiong Liu: Modified to asismilate SWE or snow depth
! 21 Jun 2019: Yeosang Yoon; Updated the file to work with the DA observation
!              space updates.
!
! !INTERFACE:
subroutine noah33_getsnowpred_PMWsnow(n, k, obs_pred)

! !USES:
  use ESMF
  use LIS_coreMod, only : LIS_rc,LIS_surface
  use noah33_lsmMod
  use PMW_snow_Mod, only : PMW_snow_struc
  use LIS_logMod,   only: LIS_logunit
  use LIS_DAobservationsMod

  implicit none
! !ARGUMENTS: 
  integer, intent(in)    :: n
  integer, intent(in)    :: k
  real                   :: obs_pred(LIS_rc%ngrid(n),LIS_rc%nensem(n))
  real                   :: snwd(LIS_rc%npatch(n,LIS_rc%lsm_index))
!EOP

  integer                :: i,t,m,gid
  real                   :: c1 !unit conversion factor

  do i=1,LIS_rc%npatch(n,LIS_rc%lsm_index),LIS_rc%nensem(n)
     do m=1,LIS_rc%nensem(n)
        t = i+m-1
        gid = LIS_surface(n,1)%tile(t)%index
        if (PMW_snow_struc(n)%dataunit .eq. 'm') then
            c1=1.0;
        elseif (PMW_snow_struc(n)%dataunit .eq. 'cm') then
            c1=100;
        elseif (PMW_snow_struc(n)%dataunit .eq. 'mm') then
            c1=1000;
        elseif (PMW_snow_struc(n)%dataunit .eq. 'inch') then
            c1=39.37;
        else
            write(LIS_logunit, *) '[ERR] unit conversion for ', &
            PMW_snow_struc(n)%dataunit, '[ERR] currently not supported'
        endif
     enddo
  enddo

  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     if (PMW_snow_struc(n)%data_var .eq. 'snow depth') then
        snwd(t) = noah33_struc(n)%noah(t)%snowh*c1   
     elseif (PMW_snow_struc(n)%data_var .eq. 'SWE') then
        snwd(t) = noah33_struc(n)%noah(t)%sneqv*c1
     else
        write(LIS_logunit, *) '[ERR] Snow assimilation of ', PMW_snow_struc(n)%data_var, &
             ' is not supported'
        return
     endif
  enddo
  
  call LIS_convertPatchSpaceToObsEnsSpace(n,k,&
       LIS_rc%lsm_index, &
       snwd,&
       obs_pred)

end subroutine noah33_getsnowpred_PMWsnow

