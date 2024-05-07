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
! !ROUTINE: noahmp36_getsnwdpred
! \label{noahmp36_getsnwdpred}
!
! !REVISION HISTORY:
! 27Feb2005: Sujay Kumar; Initial Specification
! 25Jun2006: Sujay Kumar: Updated for the ESMF design
!  02 Mar 2010: Sujay Kumar; Modified for Noah 3.1
!
! !INTERFACE:
subroutine noahmp36_getsnwdpred_TB(n, k, obs_pred)

! !USES:
  use ESMF
  use LIS_coreMod, only : LIS_rc,LIS_surface
  use noahmp36_lsmMod
  use LIS_DAobservationsMod
  use LIS_logMod,   only: LIS_logunit, LIS_verify,LIS_endrun      !kyh20170823

  implicit none
! !ARGUMENTS: 
  integer, intent(in)                                   :: n
  integer, intent(in)                                   :: k
  !real                                                 :: obs_pred(LIS_rc%ngrid(n),LIS_rc%nensem(n))  !original
  real                                                  :: obs_pred(LIS_rc%obs_ngrid(k)*4,LIS_rc%nensem(n))  !kyh20171003
  real                                                  :: snwd(LIS_rc%npatch(n,LIS_rc%lsm_index))
  real,dimension(LIS_rc%npatch(n,LIS_rc%lsm_index))     :: TB_18V, TB_18H, TB_36V, TB_36H   !kyh20171003
!EOP

  integer                :: t


!EMK...Disabled this subroutine.  The call to noahmp36_LIS_SVM is not passing
!the expected arguments, so the code is badly broken.  The task of fixing this
!rests with the original developer.
#if 1

  write(LIS_logunit,*)'[ERROR] noahmp36_getsnwdpred_TB is currently disabled.'
  write(LIS_logunit,*)'The subroutine call to noahmp36_LIS_SVM is broken'
  write(LIS_logunit,*)' (wrong arguments are passed)'
  write(LIS_logunit,*)'Edit noahmp36_getsnwdpred_TB.F90 to re-enable.'
  call LIS_endrun()

#else

  call noahmp36_LIS_SVM(n)     !kyh20170823


  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     !snwd(t) = noahmp36_struc(n)%noahmp36(t)%snowh*1000.0 !obs in mm   !original
!svk commenting out
!     TB_18V(t) = noahmp36_struc(n)%noahmp36(t)%TB_18V_lis   !K          !kyh20171003
!     TB_18H(t) = noahmp36_struc(n)%noahmp36(t)%TB_18H_lis   !K          !kyh20171003
!     TB_36V(t) = noahmp36_struc(n)%noahmp36(t)%TB_36V_lis   !K          !kyh20171003
!     TB_36H(t) = noahmp36_struc(n)%noahmp36(t)%TB_36H_lis   !K          !kyh20171003
  enddo

  !call LIS_convertPatchSpaceToObsEnsSpace(n,k,&
  !     LIS_rc%lsm_index, &
  !     snwd,&
  !     obs_pred)
  
  !---------------------------------------------kyh20171003
  call LIS_convertPatchSpaceToObsEnsSpace(n,k,&
       LIS_rc%lsm_index, &
       TB_18V,&
       obs_pred(1:LIS_rc%obs_ngrid(k),:))
  
  call LIS_convertPatchSpaceToObsEnsSpace(n,k,&
       LIS_rc%lsm_index, &
       TB_18H,&
       obs_pred(LIS_rc%obs_ngrid(k)+1:LIS_rc%obs_ngrid(k)*2,:))

  call LIS_convertPatchSpaceToObsEnsSpace(n,k,&
       LIS_rc%lsm_index, &
       TB_36V,&
       obs_pred(LIS_rc%obs_ngrid(k)*2+1:LIS_rc%obs_ngrid(k)*3,:))

  call LIS_convertPatchSpaceToObsEnsSpace(n,k,&
       LIS_rc%lsm_index, &
       TB_36H,&
       obs_pred(LIS_rc%obs_ngrid(k)*3+1:LIS_rc%obs_ngrid(k)*4,:))

  write(LIS_logunit,*) 'obs_pred=', obs_pred
  !---------------------------------------------

#endif

end subroutine noahmp36_getsnwdpred_TB

