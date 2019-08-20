!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------------
! NASA GSFC Land Data Toolkit (LDT) V1.0
!-------------------------END NOTICE -- DO NOT EDIT-----------------------------!BOP
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
