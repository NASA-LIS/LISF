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
! 
! !ROUTINE: compute_vinterp_weights
!   \label{compute_vinterp_weights} 
!
! !INTERFACE:
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
!BOP
!
! !INTERFACE: 
subroutine compute_vinterp_weights(lis_nlayers, obs_nlayers,lis_sf_d, &
  lis_rz_d, lis_d, obs_d, sf_wt, rz_wt)
      
  implicit none
! !ARGUMENTS:   
  integer   :: lis_nlayers
  integer   :: obs_nlayers
  real      :: lis_sf_d
  real      :: lis_rz_d
  real      :: lis_d(lis_nlayers)
  real      :: obs_d(obs_nlayers)
  real      :: sf_wt(obs_nlayers)
  real      :: rz_wt(obs_nlayers)
!
! !DESCRIPTION: 
!  
!  This subroutine computes vertical interpolation weights for the 
!  computation of surface and subsurface soil moisture. The routine
!  matches the soil moisture layers in the observations to the LIS 
!  model output. The interpolation weights are computed based on the
!  subset of observation layers that correspond to the model output
! 
!   The surface and root zone soil moisture (temperature) should be
!   subsequently computed as : 
!      
!     sf_sm = sum (i=1 to i=obs_nlayers) sf_wt(i)*obs_soilmoist(i)
!     rz_sm = sum (i=1 to i=obs_nlayers) rz_wt(i)*obs_soilmoist(i)
!   
!  The arguments are: 
!  \begin{description}
!   \item[lis_nlayers]      Number of layers in the LIS ouput
!   \item[obs_nlayers]      Number of layers in the observations
!   \item[lis_sf_d]         Total depth of the surface layer
!   \item[lis_rz_d]         Total depth of the subsurface layer
!   \item[lis_d]            Depths of the LIS soil layers
!   \item[obs_d]            Depths of the observation soil layers
!   \item[sf_wt]            Normalized interpolation weights for the 
!                           surface layer
!   \item[rz_wt]            Normalized interpolation weights for the 
!                           subsurface layer  
!  \end{description}
!EOP
  integer               :: i,k
  integer               :: nlayers

!computation of surface weights
  sf_wt = 0 
  rz_wt = 0 
  nlayers = 0

  do i=1, obs_nlayers
     if(obs_d(i).gt.lis_sf_d) then 
        nlayers = i 
        exit
     endif
  enddo

! compute all depths upto the last observation depth
  do k=1,nlayers
     if(obs_d(k).gt.lis_sf_d) then 
        sf_wt(k) = obs_d(k)-lis_sf_d
     else
        sf_wt(k) = obs_d(k)
     endif
  enddo
  if(sum(sf_wt).gt.0) then 
     sf_wt = sf_wt/sum(sf_wt)
  else
     sf_wt = 0
  endif

  nlayers = 0 
  do i=1, obs_nlayers
     if(obs_d(i).gt.lis_rz_d) then 
        nlayers = i 
        exit
     endif
  enddo
  
! if the deepest observation depth is shallower than the 
! LIS model depths

  if(i.gt.obs_nlayers) then 
     nlayers = i-1
  endif
! compute all depths upto the last observation depth

  do k=1,nlayers
     if(obs_d(k).gt.lis_rz_d) then 
        rz_wt(k) = obs_d(k)-lis_rz_d
     else
        rz_wt(k) = obs_d(k)
     endif
  enddo

  if(sum(rz_wt).gt.0) then 
     rz_wt = rz_wt/sum(rz_wt)
  else
     rz_wt = 0.0
  endif
  
  do i=1,obs_nlayers
     if(sf_wt(i).lt.0) then 
        print*, 'Error in surface soil weights ',sf_wt
        stop
     endif
  enddo

  do i=1,obs_nlayers
     if(rz_wt(i).lt.0) then 
        print*, 'Error in root zone soil weights ',rz_wt
        stop
     endif
  enddo
end subroutine compute_vinterp_weights
