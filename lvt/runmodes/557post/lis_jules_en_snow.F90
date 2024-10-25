!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------

! Developer note: 
! lis_jules_en_snow module is developed for creating a conceptual snowpack 
! based on LIS-JULES ensemble simulation. The steps are consistent with 
! JULES physics defined in the JULES subroutines "layersnow" and "relayersnow".
! Consistencies of snow mass and internal energy have been kept between ensemble
! means and the new snowpack.
! Shugong Wang (shugong.wang@nasa.gov) 03/10/2021 


module lis_jules_en_snow
  implicit none
  integer, parameter :: nlayer  = 3                   ! 
  integer, parameter :: nensem = 12                   !

  real, parameter, dimension(3) :: dzsnow= [0.04, 0.12, 0.34]

  real, parameter :: hcapi               = 2100.0
  real, parameter :: hcapw               = 4180.0
  real, parameter :: tm                  = 273.15
  
  integer, parameter ::     i_relayer_opt = 1          ! i_relayer_opt is either 0 (ip_relayer_linear) or 1 (ip_relayer_rgrain_inv)
                                                       ! i_relayer_opt should be consistent to JULES snow name list file
  integer, parameter ::     ip_relayer_linear = 0      ! *** relayer based on thickness  
  integer, parameter ::     ip_relayer_rgrain_inv = 1  ! *** as the linear scheme, but relayering the inverse of the grain size
  
  ! input part, populate the following input variables before calling the step subroutines 
  real, dimension(nensem, nlayer) :: en_sice           
  real, dimension(nensem, nlayer) :: en_sliq
  real, dimension(nensem, nlayer) :: en_tsnow
  real, dimension(nensem, nlayer) :: en_rgrainl
  real, dimension(nensem, nlayer) :: en_ds 
  integer, dimension(nensem)      :: en_nsnow


  ! local part 
  real, dimension(nensem, nlayer) :: en_C ! heat capacity (J/kg/K)
  real, dimension(nensem, nlayer) :: en_e ! internal energy (J/m2)

  real, dimension(0:nlayer) :: sice0      ! average ice water content across ensemble for each snow laye
  real, dimension(0:nlayer) :: sliq0      ! average liquid water content across ensemble for each snow layer
  real, dimension(0:nlayer) :: e0         ! average internal energy across ensemble for each snow layer
  real, dimension(0:nlayer) :: r0         ! weighted average of grain size across ensemble for each snow layer
  real, dimension(0:nlayer) :: d0         ! average depth across ensemble for each snow layer
  real, dimension(nlayer)   :: u0         ! layer energy contents.
  real, dimension(nlayer)   :: newremains ! available (unfilled) depth in new layer (m).
  real                      :: D          ! average snow depth: D = sum of d0 over snow layers
  real, parameter :: thin_snow_limit = 1.0e-6

! output part, grab values in these variables after calling the step subroutines 
  real, dimension(nlayer) :: new_ds
  real, dimension(nlayer) :: new_sice
  real, dimension(nlayer) :: new_sliq
  real, dimension(nlayer) :: new_tsnow
  real, dimension(nlayer) :: new_rgrainl
  real, dimension(nlayer) :: new_rho_snow
  real                    :: new_snowmass
  real                    :: new_rgrain
  real                    :: new_rho_snow_grnd
  real                    :: new_snowdepth
  real                    :: new_rho_grnd
  integer                 :: new_nsnow 

contains
  real function layer_mean(en_var_l, l)
    implicit none
    real, dimension(nensem, nlayer) ,intent(in):: en_var_l
    integer :: l
    ! local varialbes 
    integer :: n
    real :: s
    
    s = 0.0
    do n =1, nensem
      s = s + en_var_l(n, l)
    end do
    layer_mean = s/nensem
    
  end function

  !1. Calculate average ice water content across ensemble for each snow layer: sice0 = 1/N sum(sice)
  subroutine step_1()
    implicit none 
    integer :: l
    
    sice0(:) = 0.0
    do l=1, nlayer
      sice0(l) = layer_mean(en_sice, l)
    end do
    
  end subroutine step_1
  
  !2. Calculate average liquid water content across ensemble for each snow layer: sliq0 = 1/N sum(sliq)
  subroutine step_2()
    implicit none 
    integer :: l
    
    sliq0(:) = 0.0
    do l=1, nlayer
      sliq0(l) = layer_mean(en_sliq, l)
    end do

  end subroutine step_2
  
  !3. For each snow layer: calculate heat capacity for each ensemble member: C = sice * Ci + sliq * Cw
  subroutine step_3()
    implicit none
    integer :: l, n

    en_C(:,:) = 0.0
    do n=1, nensem
      do l=1, nlayer
        en_C(n, l) = en_sice(n, l)*hcapi + en_sliq(n, l)*hcapw
      end do
    end do

  end subroutine step_3

  !4. For each snow layer: calculate the energy for each ensemble member: e = C(T-Tw)
  !*** Tw is assumed to be tm (temperature at which fresh water freezes and ice melts)
  subroutine step_4()
    implicit none
    integer :: l, n

    en_e(:,:) = 0.0
    do n=1, nensem
      do l=1, nlayer
        en_e(n, l) = en_C(n, l)*(en_tsnow(n, l) - tm)
      end do
    end do

  end subroutine step_4

  !5. Calculate average energy across ensemble for each snow layer: e0 = 1/N sum(e)
  subroutine step_5()
    implicit none 
    integer :: l
    
    e0(:) = 0.0 
    do l=1, nlayer
      e0(l) = layer_mean(en_e, l)
    end do

  end subroutine step_5

  !6. Calculate weighted average of grain size across ensemble for each snow layer:  
  !   r0 = 1/N sum( (sice+sliq)*r / (sice0+sliq0)) 
  subroutine step_6()
    implicit none
    integer :: l, n, n_empty
   
    r0(:) = 0.0
    do l=1, nlayer
      n_empty = 0
      do n=1, nensem
        if (sice0(l)+sliq0(0) > 0.0) then
          r0(l) = r0(l) + (en_sice(n,l)+en_sliq(n,l))*en_rgrainl(n,l)/(sice0(l)+sliq0(0))
        else
          n_empty = n_empty + 1
        endif
      end do
      if (nensem > n_empty) then
        r0(l) = r0(l)/(nensem - n_empty)
      else
        r0(l) = 50.0 ! 
      endif
    end do
  end subroutine step_6

  !7. Calculate average depth across ensemble for each snow layer: d0 = 1/N sum(dze)
  !8. Calculate average snow depth: D = sum of d0 over snow layers
  subroutine step_7_8()
    implicit none 
    integer :: l

    D = 0.0
    do l=1, nlayer
      d0(l) = layer_mean(en_ds, l)
      D = D + d0(l)
    end do
  end subroutine step_7_8

  !9. Calculate new snow layers depths using D as input into JULES routine layersnow (outputs new snow layer thicknesses: dz)
  subroutine step_9()
    implicit none
    integer :: l
    real :: remains
    
    new_ds(:) = 0.0
    ! only divide snowpack into layers if depth is >= a threshold
    if (D >= dzsnow(1)) then
      remains = D
      do l = 1, nlayer
        new_ds(l) = dzsnow(l)
        
        ! set number of snow layers
        new_nsnow = l 
        
        remains = remains - dzsnow(l)
        if (remains <= dzsnow(l) .or. l==nlayer) then
          new_ds(l) = new_ds(l) + remains
          exit
        end if

      end do
    else
      if (D>0) then
        new_nsnow=1
        new_ds(1) = D
      endif
    endif

    new_snowdepth = sum(new_ds)

  end subroutine step_9

  !10. Calculate new snow layer properties using the method in the JULES routine relayersnow from lines: 323 - 440
  subroutine step_10
    implicit none
    real    :: csnow
    real    :: oldremains ! remaining depth in an old layer (m).
    real    :: wt         ! weight given to a layer value.

    integer :: l, nold, new, old, iznew, izz 
    real, dimension(nlayer) :: u
    
    ! number of (effective) layers before adjustment. 
    nold = maxval(en_nsnow) 


    !!! initialize accumulations for new layer values  
    u(:)            = 0.0
    new_sice(:)     = 0.0
    new_sliq(:)     = 0.0
    new_rgrainl(:)  = 0.0
    
    !!! initialize with all new layers empty 
    do l = 1, new_nsnow
      newremains(l) = new_ds(l) ! ! snow layer thicknesses (m), new
    end do


    iznew = 1 
    ! loop over the old layers
    ! 0 represent new snow. we set an empty layer 0 for this case
    do old=0, nold
      ! all of this old layer remains to be reassigned to new layer(s).
      oldremains = d0(old)

      ! point to first new layer with remaining space.
      izz = iznew 
      ! loop over new layers with remaining space
      do new = izz, new_nsnow
        if (oldremains>newremains(new)) then ! the old remains is more than the capacity of the current new layer
          ! The remaining depth in the new layer will be exhausted by some or
          ! all of the remaining depth from the old layer.
          ! Note: newremains <-> left capacity/depth of the new layer
          !       oldremains <-> remain snow of the old layer
          ! decrement old layer by the remaining space in new layer.
          oldremains = oldremains - newremains(new)

          ! add properties from old layers to accumulation for new layer 
          ! note that wt is <=1 since here we have oldremains > newremains
          ! and oldremains<=d0 ???
          if ( d0(old) > thin_snow_limit) then
            wt        = newremains(new)/d0(old)
            u(new)    = u(new) + e0(old)*wt
            new_sice(new) = new_sice(new) + sice0(old) * wt
            new_sliq(new) = new_sliq(new) + sliq0(old) * wt
            
            select case(i_relayer_opt)
            case (ip_relayer_linear)
              new_rgrainl(new) = new_rgrainl(new) + r0(old)*newremains(new)
            case (ip_relayer_rgrain_inv)
              new_rgrainl(new) = new_rgrainl(new) + newremains(new)/r0(old)
            end select
          endif
        
          ! update the pointer to the next new layer with space
          izz = new + 1

        else !  the old layer will be exhausted by this increment.
          ! decrement available space in the new layer.
          newremains(new) = newremains(new) - oldremains

          ! add properties from old layer to accumulation for new layer
          if (d0(old) > thin_snow_limit) then 
            wt     = oldremains/d0(old)
            u(new) = u(new) + e0(old)*wt
            new_sice(new) = new_sice(new) + sice0(old) * wt
            new_sliq(new) = new_sliq(new) + sliq0(old) * wt
            select case(i_relayer_opt)
            case (ip_relayer_linear)
              new_rgrainl(new) = new_rgrainl(new) + r0(old)*oldremains
            case (ip_relayer_rgrain_inv)
              new_rgrainl(new) = new_rgrainl(new) + oldremains/r0(old)
            end select
          endif  
          ! proceed to the next old layer by exiting from the new layer loop
          exit 
        end if
      end do ! new layers 

      ! update pointer to the next layer with space
      iznew = izz 
    end do

    ! diagnose layer temperatures and densities 
    do l =1, new_nsnow
      csnow = new_sice(l) * hcapi + new_sliq(l) * hcapw
      new_tsnow(l) = tm + u(l)/csnow
      new_rho_snow(l) = (new_sice(l) + new_sliq(l))/new_ds(l)
      select case(i_relayer_opt)
      case (ip_relayer_linear)
        new_rgrainl(l) = new_rgrainl(l)/new_ds(l)
      case (ip_relayer_rgrain_inv)
        new_rgrainl(l) = new_ds(l)/new_rgrainl(l)
      end select
    end do

    ! snow surface grain size for radiative calculations
    new_rgrain = new_rgrainl(1)
    

  end subroutine step_10

  !11. Calculate snowmass (sum of ice and liquid water contents)
  subroutine step_11
    implicit none
    integer :: l 

    new_snowmass = 0.0
    do l =1, new_nsnow
      new_snowmass = new_snowmass + new_sice(l) + new_sliq(l)
    end do

    ! diagnose bulk density of snowpack 
    new_rho_grnd = new_snowmass/new_snowdepth
  end subroutine step_11
end module lis_jules_en_snow

