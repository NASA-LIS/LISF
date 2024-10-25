!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
! 
! !USES:   
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
subroutine vinterp(nlis, nobs, udef, lis_data, obs_data, wts)

  implicit none

  integer        :: nlis
  integer        :: nobs
  real           :: udef
  real           :: lis_data(nlis)
  real           :: obs_data(nobs)
  real           :: wts(nlis,nobs)

  integer        :: k,j

  lis_data = 0 

  do k=1,nlis
     do j=1,nobs
        if(obs_data(k).ne.udef) then 
           lis_data(k) = lis_data(k) + obs_data(k)*wts(k,j)
        endif
     enddo
  enddo

  
end subroutine vinterp
