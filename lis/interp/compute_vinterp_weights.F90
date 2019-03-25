!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------------
! NASA GSFC Land surface Verification Toolkit (LVT) V1.0
!-------------------------END NOTICE -- DO NOT EDIT-----------------------------
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
subroutine compute_vinterp_weights(obs_nlayers,lis_sf_d, &
  lis_rz_d, obs_d, sf_wt, rz_wt)
      
  implicit none
! !ARGUMENTS:   
  integer   :: obs_nlayers
  real      :: lis_sf_d
  real      :: lis_rz_d
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
!  
!  Algorithm 2 (from NASA MSFC): 
! 
!  Examine mid-points between successive obs levels. 
!  Assuming standard SCAN depths of 5, 10, 20, 50, and 100cm, 
!  the mid points are 7.5, 15, 35, and 75 cm. 
!   Top level has special handling for relative weight assignment, 
!   simply being the average of top two obs levels (i.e., 7.5cm)
!   Intermediate level weights for the handling of 0-1m root zone 
!   are half the distance of the midpoint below and above a given 
!   obs level.  E.G. the 20cm obs level will be assigned a weight 
!   of 0.5*(35-15), or 10.0. 
!   If an obs level is BELOW the bottom of the surface or root zone, 
!   examine the depth of the midpoint above and below the bottom zone 
!   layer.  If the mid point is below the bottom of the zone, assign 
!   0 weight.  Otherwise, assign a relative weight of zone depth 
!   minus last mid point above base of root zone.  E.G. For 10cm 
!   surface in model and standard SCAN obs levels, bottom 
!   relative weight is 10 – 0.5*(10+5)) = 10 - 7.5 = 2.5. 
! 
!  The arguments are: 
!  \begin{description}
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
  
  nlayers = 1

  do i=1, obs_nlayers
     if(obs_d(i).gt.lis_sf_d) then 
        nlayers = i 
        exit
     endif
  enddo
  
  nlayers = max(1,nlayers-1)
  ! compute all depths upto the last observation depth
  
  ! first check if the top layer is too thick (difference > 2cm) then 
  ! do not use for surface soil moisture
  if(obs_d(1)-lis_sf_d > 0.02) then 
     sf_wt = 0 
  else
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
  endif
  
  nlayers = 1
  do i=1, obs_nlayers
     if(obs_d(i).gt.lis_rz_d) then 
        nlayers = i 
        exit
     endif
  enddo
  nlayers = max(1,nlayers-1)
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
