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
! !ROUTINE: clsmf25_gettwspred
! \label{clsmf25_gettwspred}
!
! !REVISION HISTORY:
! 27Feb2005: Sujay Kumar; Initial Specification
! 25Jun2006: Sujay Kumar: Updated for the ESMF design
! 28Sep2011: Ben Zaitchik: Applied to GRACE
!
! !INTERFACE:
subroutine clsmf25_gettwspred(n,k,obs_pred)
! !USES:
  use ESMF
  use LIS_coreMod
  use clsmf25_lsmMod
  use LIS_logMod
  use clsmf25_tws_DAlogMod
  use LIS_DAobservationsMod
        
!EOP

  implicit none
! !ARGUMENTS: 
  integer, intent(in)    :: n
  integer, intent(in)    :: k
!  real                   :: clmnwater(3,LIS_rc%obs_ngrid(k),LIS_rc%nensem(n))
  real                   :: obs_pred(LIS_rc%obs_ngrid(k),LIS_rc%nensem(n))
!
! !DESCRIPTION:
!
!  Returns the TWS obs pred (model's estimate of 
!  observations) for data assimilation
! 
!  The arguments are: 
!  \begin{description}
!  \item[n] index of the nest \newline
!  \item[obs\_pred] model's estimate of observations \newline
!  \end{description}
!EOP
  integer                :: t
  real                   :: tws(LIS_rc%npatch(n,LIS_rc%lsm_index))

 do t=1, LIS_rc%npatch(n,LIS_rc%lsm_index)
     tws(t) = (CLSMpred_struc(n)%clmnwater(1,t) + &	
	 CLSMpred_struc(n)%clmnwater(2,t)+&
         CLSMpred_struc(n)%clmnwater(3,t))/3.	
  enddo
  call LIS_convertPatchSpaceToObsEnsSpace(n,k,&
       LIS_rc%lsm_index, &
       tws,&
       obs_pred)

end subroutine clsmf25_gettwspred

