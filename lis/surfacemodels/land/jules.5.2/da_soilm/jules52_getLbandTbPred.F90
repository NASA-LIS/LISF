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
! !ROUTINE: jules52_getLbandTbPred
! \label{jules52_getLbandTbPred}
!
! !REVISION HISTORY:
! 13 Sept 2012: Sujay Kumar; Initial Specification
! 20 Dec 2018: Mahdi Navari; Modified for JULES 5.2
!
! !INTERFACE:
subroutine jules52_getLbandTbPred(n, obs_pred)
! !USES:
  use ESMF
  use LIS_coreMod, only : LIS_rc,LIS_surface
  use LIS_RTMMod,  only : LIS_forwardState, LIS_RTM_run
  use LIS_logMod,  only : LIS_verify
  
  use jules52_lsmMod
!EOP

  implicit none
! !ARGUMENTS: 
  integer, intent(in)    :: n
  real                   :: obs_pred(LIS_rc%ngrid(n)*2,LIS_rc%nensem(n))
  integer                :: count1(LIS_rc%ngrid(n)*2,LIS_rc%nensem(n))
!
! !DESCRIPTION:
!
!  Returns the Soil moisture obs pred (model's estimate of 
!  observations) for data assimilation
! 
!  The arguments are: 
!  \begin{description}
!  \item[n] index of the nest \newline
!  \item[obs\_pred] model's estimate of observations \newline
!  \end{description}
!EOP
  integer                :: i,t,m,gid
  real, pointer          :: TbH(:), TbV(:)
  type(ESMF_Field)       :: varField
  integer                :: status

!call the forward model
  call LIS_RTM_run(n)

  call ESMF_StateGet(LIS_forwardState(n), "Lband_Hpol", varField, rc=status)
  call LIS_verify(status, &
       "Error in StateGet in jules52_getLbandTbPred for Lband_Hpol")

  call ESMF_FieldGet(varField, localDE=0,farrayPtr=TbH, rc=status)
  call LIS_verify(status, &
       'Error in FieldGet in jules52_getLbandTbPred for Lband_Hpol')

  call ESMF_StateGet(LIS_forwardState(n), "Lband_Vpol", varField, rc=status)
  call LIS_verify(status, &
       "Error in StateGet in jules52_getLbandTbPred for Lband_Vpol")

  call ESMF_FieldGet(varField, localDE=0,farrayPtr=TbV, rc=status)
  call LIS_verify(status, &
       'Error in FieldGet in jules52_getLbandTbPred for Lband_Vpol')
 
  obs_pred = 0.0
  count1 = 0.0 
  do i=1,LIS_rc%npatch(n,LIS_rc%lsm_index),LIS_rc%nensem(n)
     do m=1,LIS_rc%nensem(n)
        t = i+m-1
        gid = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%index        
        obs_pred(gid,m)= obs_pred(gid,m) + TbH(t)
        count1(gid,m) = count1(gid,m) +1
     enddo
  enddo

  do i=1,LIS_rc%npatch(n,LIS_rc%lsm_index),LIS_rc%nensem(n)
     do m=1,LIS_rc%nensem(n)
        t = i+m-1
        gid = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%index + &
             LIS_rc%ngrid(n)
        obs_pred(gid,m)= obs_pred(gid,m) + TbV(t)
        count1(gid,m) = count1(gid,m) +1
     enddo
  enddo
  
  do i=1,LIS_rc%ngrid(n)*2
     do m=1,LIS_rc%nensem(n)
        obs_pred(i,m) = obs_pred(i,m)/(count1(i,m))
     enddo
  enddo

end subroutine jules52_getLbandTbPred

