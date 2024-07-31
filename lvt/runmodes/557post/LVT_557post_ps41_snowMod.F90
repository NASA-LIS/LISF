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
! MODULE: LVT_557post_ps41_snowMod
!
! PURPOSE:
!
! This module was developed to create a conceptual snowpack based on a
! LIS-JULES ensemble simulation with "PS41" physics, to pass to the GALWEM
! NWP model which is coupled with JULES. The steps are consistent
! with JULES physics defined in the JULES subroutines "layersnow" and
! "relayersnow".  Consistencies of snow mass and internal energy have been
! kept between ensemble means and the new snowpack.
!
! REVISION HISTORY:
! 10 Mar 2021: Shugong Wang, Initial Specification.
! 19 Mar 2021: Eric Kemp, Ported to LVT.
!------------------------------------------------------------------------------

module LVT_557post_ps41_snowMod

  ! Defaults
  implicit none
  private

  integer, parameter :: nlayer  = 3 ! JULES PS41 uses three-layer snow physics

  real, parameter, dimension(3) :: dzsnow= (/ 0.04, 0.12, 0.34 /)

  real, parameter :: hcapi               = 2100.0
  real, parameter :: hcapw               = 4180.0
  real, parameter :: tm                  = 273.15

  ! i_relayer_opt is either 0 (ip_relayer_linear) or 1 (ip_relayer_rgrain_inv)
  ! i_relayer_opt should be consistent to JULES snow name list file
  integer, parameter :: i_relayer_opt = 1

  ! *** relayer based on thickness
  integer, parameter :: ip_relayer_linear = 0

  ! *** as the linear scheme, but relayering the inverse of the grain size
  integer, parameter :: ip_relayer_rgrain_inv = 1

  ! Input part, populate the following input variables before calling the
  ! step subroutines
  real, allocatable :: en_sice(:,:)
  real, allocatable :: en_sliq(:,:)
  real, allocatable :: en_tsnow(:,:)
  real, allocatable :: en_rgrainl(:,:)
  real, allocatable :: en_ds(:,:)
  integer, allocatable :: en_nsnow(:)

  ! Output part, grab values in these variables after calling the step
  ! subroutines
  real, dimension(nlayer) :: new_ds
  real, dimension(nlayer) :: new_sice
  real, dimension(nlayer) :: new_sliq
  real, dimension(nlayer) :: new_tsnow
  real, dimension(nlayer) :: new_rgrainl
  real, dimension(nlayer) :: new_rho_snow
  real :: new_snowmass
  real :: new_rgrain
  real :: new_snowdepth
  real :: new_rho_grnd
  integer :: new_nsnow

  ! Internal variables
  real, allocatable :: en_C(:,:) ! heat capacity (J/kg/K)
  real, allocatable :: en_e(:,:) ! internal energy (J/m2)

  ! Average ice water content across ensemble for each snow layer
  real, dimension(0:nlayer) :: sice0
  ! Average liquid water content across ensemble for each snow layer
  real, dimension(0:nlayer) :: sliq0
  ! Average internal energy across ensemble for each snow layer
  real, dimension(0:nlayer) :: e0
  ! Weighted average of grain size across ensemble for each snow layer
  real, dimension(0:nlayer) :: r0
  ! Average depth across ensemble for each snow layer
  real, dimension(0:nlayer) :: d0
  ! Layer energy contents.
  real, dimension(nlayer)   :: u0
  ! Available (unfilled) depth in new layer (m).
  real, dimension(nlayer)   :: newremains
  ! Average snow depth: D = sum of d0 over snow layers
  real                      :: D

  real, parameter :: thin_snow_limit = 1.0e-6

  ! Internal pointers for referencing PS41 JULES fields in LVT dataEntry
  ! linked list.
  real, pointer, save :: snowIce(:,:,:)
  real, pointer, save :: snowLiq(:,:,:)
  real, pointer, save :: snowTProf(:,:,:)
  real, pointer, save :: layerSnowGrain(:,:,:)
  real, pointer, save :: layerSnowDepth(:,:,:)
  real, pointer, save :: ActSnowNL(:,:,:)

  ! Internal arrays for storing processed multi-layer snow
  real, allocatable, save :: actSnowNL_final(:)
  real, allocatable, save :: grndSnow_final(:)
  real, allocatable, save :: layerSnowDensity_final(:,:)
  real, allocatable, save :: layerSnowDepth_final(:,:)
  real, allocatable, save :: layerSnowGrain_final(:,:)
  real, allocatable, save :: SWE_final(:)
  real, allocatable, save :: snowDensity_final(:)
  real, allocatable, save :: snowDepth_final(:)
  real, allocatable, save :: snowGrain_final(:)
  real, allocatable, save :: snowIce_final(:,:)
  real, allocatable, save :: snowLiq_final(:,:)
  real, allocatable, save :: snowTProf_final(:,:)
  real, allocatable, save :: surftSnow_final(:)

  ! Public routines
  public :: LVT_init_jules_ps41_ens_snow
  public :: LVT_prep_jules_ps41_ens_snow_var
  public :: LVT_proc_jules_ps41_ens_snow
  public :: LVT_fetch_jules_ps41_ens_snow_final
  public :: LVT_cleanup_jules_ps41_ens_snow
  public :: LVT_set_SWE_metadata
  public :: LVT_set_SnowDensity_metadata
  public :: LVT_set_LayerSnowDensity_metadata
  public :: LVT_set_SnowGrain_metadata
  public :: LVT_set_SnowDepth_metadata
  public :: LVT_set_GrndSnow_metadata
  public :: LVT_set_SurftSnow_metadata
  public :: LVT_deallocate_metadata

contains

  ! Allocate all step arrays based on number of ensemble members.
  ! Note that nlayers is hardwired since this targets PS41.
  subroutine allocate_step_arrays(nensem)
    implicit none
    integer, intent(in) :: nensem
    allocate(en_sice(nensem, nlayer))
    allocate(en_sliq(nensem, nlayer))
    allocate(en_tsnow(nensem, nlayer))
    allocate(en_rgrainl(nensem, nlayer))
    allocate(en_ds(nensem, nlayer))
    allocate(en_nsnow(nensem))
    allocate(en_C(nensem, nlayer))
    allocate(en_e(nensem, nlayer))
  end subroutine allocate_step_arrays

  ! Allocate arrays for final data
  subroutine allocate_final_arrays(nc, nr)
    implicit none
    integer, intent(in) :: nc
    integer, intent(in) :: nr
    allocate(actSnowNL_final(nc*nr))
    allocate(grndSnow_final(nc*nr))
    allocate(layerSnowDensity_final(nc*nr, nlayer))
    allocate(layerSnowDepth_final(nc*nr, nlayer))
    allocate(layerSnowGrain_final(nc*nr, nlayer))
    allocate(SWE_final(nc*nr))
    allocate(snowDensity_final(nc*nr))
    allocate(snowDepth_final(nc*nr))
    allocate(snowGrain_final(nc*nr))
    allocate(snowIce_final(nc*nr, nlayer))
    allocate(snowLiq_final(nc*nr, nlayer))
    allocate(snowTProf_final(nc*nr, nlayer))
    allocate(surftSnow_final(nc*nr))
  end subroutine allocate_final_arrays

  ! Initialize routine for JULES PS41 ensemble snow
  subroutine LVT_init_jules_ps41_ens_snow()
    implicit none
    call nullify_pointers()
  end subroutine LVT_init_jules_ps41_ens_snow

  ! Nullify all internal pointers
  subroutine nullify_pointers()
    implicit none
    nullify(snowIce)
    nullify(snowLiq)
    nullify(snowTProf)
    nullify(layerSnowGrain)
    nullify(layerSnowDepth)
    nullify(ActSnowNL)
  end subroutine nullify_pointers

  ! Deallocate memory for step routines.  Final arrays are preserved for
  ! later extraction by LVT.
  subroutine deallocate_step_arrays()
    implicit none
    if (allocated(en_sice)) deallocate(en_sice)
    if (allocated(en_sliq)) deallocate(en_sliq)
    if (allocated(en_tsnow)) deallocate(en_tsnow)
    if (allocated(en_rgrainl)) deallocate(en_rgrainl)
    if (allocated(en_ds)) deallocate(en_ds)
    if (allocated(en_nsnow)) deallocate(en_nsnow)
    if (allocated(en_C)) deallocate(en_C)
    if (allocated(en_e)) deallocate(en_e)
  end subroutine deallocate_step_arrays

  ! Calculate mean snow layer depth from ensemble members
  real function layer_mean(nensem, en_var_l, l)
    implicit none
    integer, intent(in) :: nensem
    real, dimension(nensem, nlayer) ,intent(in):: en_var_l
    integer, intent(in) :: l
    integer :: n
    real :: s
    s = 0.0
    do n = 1, nensem
      s = s + en_var_l(n,l)
    end do
    layer_mean = s / nensem
  end function layer_mean

  !1. Calculate average ice water content across ensemble for each snow
  !   layer: sice0 = 1/N sum(sice)
  subroutine step_1(nensem)
    implicit none
    integer, intent(in) :: nensem
    integer :: l
    sice0(:) = 0.0
    do l = 1, nlayer
      sice0(l) = layer_mean(nensem, en_sice, l)
    end do
  end subroutine step_1

  !2. Calculate average liquid water content across ensemble for each snow
  !   layer: sliq0 = 1/N sum(sliq)
  subroutine step_2(nensem)
    implicit none
    integer, intent(in) :: nensem
    integer :: l
    sliq0(:) = 0.0
    do l = 1, nlayer
      sliq0(l) = layer_mean(nensem, en_sliq, l)
    end do
  end subroutine step_2

  !3. For each snow layer: calculate heat capacity for each ensemble member:
  !   C = sice * Ci + sliq * Cw
  subroutine step_3(nensem)
    implicit none
    integer, intent(in) :: nensem
    integer :: l, n
    en_C(:,:) = 0.0
    do l = 1, nlayer
       do n = 1, nensem
         en_C(n,l) = en_sice(n,l)*hcapi + en_sliq(n,l)*hcapw
      end do
    end do
  end subroutine step_3

  !4. For each snow layer: calculate the energy for each ensemble member:
  !   e = C(T-Tw)
  !   *** Tw is assumed to be tm (temperature at which fresh water freezes and
  !   ice melts)
  subroutine step_4(nensem)
    implicit none
    integer, intent(in) :: nensem
    integer :: l, n
    en_e(:,:) = 0.0
    do l = 1, nlayer
       do n = 1, nensem
         en_e(n,l) = en_C(n,l)*(en_tsnow(n,l) - tm)
      end do
    end do
  end subroutine step_4

  !5. Calculate average energy across ensemble for each snow layer:
  !   e0 = 1/N sum(e)
  subroutine step_5(nensem)
    implicit none
    integer, intent(in) :: nensem
    integer :: l
    e0(:) = 0.0
    do l = 1, nlayer
      e0(l) = layer_mean(nensem, en_e, l)
    end do
  end subroutine step_5

  !6. Calculate weighted average of grain size across ensemble for each snow
  !   layer:  r0 = 1/N sum( (sice+sliq)*r / (sice0+sliq0))
  subroutine step_6(nensem)
    implicit none
    integer, intent(in) :: nensem
    integer :: l, n, n_empty
    r0(:) = 0.0
    do l = 1, nlayer
      n_empty = 0
      do n = 1, nensem
        if (sice0(l)+sliq0(l) > 0.0) then
           r0(l) = r0(l) + &
                (en_sice(n,l) + en_sliq(n,l)) * &
                en_rgrainl(n,l) / (sice0(l) + sliq0(l))
        else
          n_empty = n_empty + 1
        endif
      end do
      if (nensem > n_empty) then
        r0(l) = r0(l) / (nensem - n_empty)
      else
        r0(l) = 50.0 !
      endif
    end do
  end subroutine step_6

  !7. Calculate average depth across ensemble for each snow layer:
  !   d0 = 1/N sum(dze)
  !8. Calculate average snow depth: D = sum of d0 over snow layers
  subroutine step_7_8(nensem)
    implicit none
    integer, intent(in) :: nensem
    integer :: l
    D = 0.0
    do l = 1, nlayer
      d0(l) = layer_mean(nensem, en_ds, l)
      D = D + d0(l)
    end do
  end subroutine step_7_8

  !9. Calculate new snow layers depths using D as input into JULES routine
  !   layersnow (outputs new snow layer thicknesses: dz)
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
        if (remains <= dzsnow(l) .or. l == nlayer) then
          new_ds(l) = new_ds(l) + remains
          exit
        end if
      end do
    else
      if (D > 0) then
        new_nsnow=1
        new_ds(1) = D
      endif
    endif

    new_snowdepth = sum(new_ds)

  end subroutine step_9

  !10. Calculate new snow layer properties using the method in the JULES
  !    routine relayersnow from lines: 323 - 440
  subroutine step_10()

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

    ! EMK...Other initializations
    new_rho_snow(:) = 0.0
    new_tsnow(:) = tm

    !!! initialize with all new layers empty
    do l = 1, new_nsnow
      newremains(l) = new_ds(l) ! ! snow layer thicknesses (m), new
    end do

    iznew = 1
    ! loop over the old layers
    ! 0 represent new snow. we set an empty layer 0 for this case
    do old = 0, nold
      ! all of this old layer remains to be reassigned to new layer(s).
      oldremains = d0(old)

      ! point to first new layer with remaining space.
      izz = iznew
      ! loop over new layers with remaining space
      do new = izz, new_nsnow
         if (oldremains  > newremains(new)) then
          ! the old remains is more than the capacity of the current new layer
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
            wt        = newremains(new) / d0(old)
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
            wt     = oldremains / d0(old)
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
    do l = 1, new_nsnow
      csnow = new_sice(l) * hcapi + new_sliq(l) * hcapw
      new_tsnow(l) = tm + u(l)/csnow
      new_rho_snow(l) = (new_sice(l) + new_sliq(l)) / new_ds(l)
      !EMK...Add limits to snow density
      new_rho_snow(l) = max(109., min(917., new_rho_snow(l)))

      select case(i_relayer_opt)
      case (ip_relayer_linear)
        new_rgrainl(l) = new_rgrainl(l) / new_ds(l)
      case (ip_relayer_rgrain_inv)
        new_rgrainl(l) = new_ds(l) / new_rgrainl(l)
      end select
    end do

    ! EMK...Enforce snow grain size limits
    do l = 1, nlayer
       new_rgrainl(l) = min(2000., max(new_rgrainl(l), 50.) )
    end do

    ! snow surface grain size for radiative calculations
    new_rgrain = new_rgrainl(1)

  end subroutine step_10

  !11. Calculate snowmass (sum of ice and liquid water contents)
  subroutine step_11()
    implicit none
    integer :: l
    new_snowmass = 0.0
    do l = 1, new_nsnow
      new_snowmass = new_snowmass + new_sice(l) + new_sliq(l)
    end do
    ! diagnose bulk density of snowpack
    new_rho_grnd = new_snowmass / new_snowdepth
    ! EMK...Add limits to bulk density
    new_rho_grnd = max(109., min(917., new_rho_grnd))
  end subroutine step_11

  ! Set pointer to passed array based on variable name
  subroutine LVT_prep_jules_ps41_ens_snow_var(short_name, vlevels, data, &
       is_ps41_snow_var)

    ! Defaults
    implicit none

    ! Arguments
    character(len=*), intent(in) :: short_name
    integer, intent(in) :: vlevels
    real, target, intent(in) :: data(:,:,:)
    logical, intent(out) :: is_ps41_snow_var

    is_ps41_snow_var = .true. ! First guess

    ! Preliminary check for 3 vertical levels (used with PS41 snow fields).
    ! This doesn't apply for ActSnowNL, since that is a single layer variable.
    select case(trim(short_name))
    case ("ActSnowNL_inst")
       continue
    case default
       if (vlevels .ne. 3) then
          is_ps41_snow_var = .false.
          return
       end if
    end select

    ! Now update the appropriate pointer
    select case(trim(short_name))
    case ("SnowIce_inst")
       snowIce => data
    case ("SnowLiq_inst")
       snowLiq => data
    case ("SnowTProf_inst")
       snowTProf => data
    case ("LayerSnowGrain_inst")
       layerSnowGrain => data
    case ("LayerSnowDepth_inst")
       layerSnowDepth => data
    case ("ActSnowNL_inst")
       actSnowNL => data
    case default
       is_ps41_snow_var = .false.
    end select

  end subroutine LVT_prep_jules_ps41_ens_snow_var

  ! Process the JULES PS41 snow variables
  subroutine LVT_proc_jules_ps41_ens_snow()

    ! Modules
    use LVT_coreMod, only: LVT_domain, LVT_rc

    ! Defaults
    implicit none

    ! Locals
    integer :: c, r, m, k, gid, ij

    ! Initializations
    call allocate_step_arrays(LVT_rc%nensem)
    call allocate_final_arrays(LVT_rc%lnc, LVT_rc%lnr)
    actSnowNL_final = LVT_rc%udef
    grndSnow_final = LVT_rc%udef
    layerSnowDensity_final = LVT_rc%udef
    layerSnowDepth_final = LVT_rc%udef
    layerSnowGrain_final = LVT_rc%udef
    SWE_final = LVT_rc%udef
    snowDensity_final = LVT_rc%udef
    snowDepth_final = LVT_rc%udef
    snowGrain_final = LVT_rc%udef
    snowIce_final = LVT_rc%udef
    snowLiq_final = LVT_rc%udef
    snowTProf_final = LVT_rc%udef
    surftSnow_final = LVT_rc%udef

    ! For each land point, apply PS41 ensemble post-processing
    do r = 1, LVT_rc%lnr
       do c = 1, LVT_rc%lnc
          gid = LVT_domain%gindex(c,r)
          if (gid == -1) cycle

          ij = c + (r-1)*LVT_rc%lnc

          ! Load ensemble members into step arrays
          do k = 1, nlayer
             do m = 1, LVT_rc%nensem
                en_sice(m,k) = snowIce(gid,m,k)
                en_sliq(m,k) = snowLiq(gid,m,k)
                en_tsnow(m,k) = snowTProf(gid,m,k)
                en_rgrainl(m,k) = layerSnowGrain(gid,m,k)
                en_ds(m,k) = layerSnowDepth(gid,m,k)
             end do ! m
          end do ! k
          do m = 1, LVT_rc%nensem
             en_nsnow(m) = actSnowNL(gid,m,1)
          end do ! m

          ! Execute most of the relayering algorithm
          call step_1(LVT_rc%nensem)
          call step_2(LVT_rc%nensem)
          call step_3(LVT_rc%nensem)
          call step_4(LVT_rc%nensem)
          call step_5(LVT_rc%nensem)
          call step_6(LVT_rc%nensem)
          call step_7_8(LVT_rc%nensem)
          call step_9()

          ! Handle no-snow case
          if (new_snowdepth .eq. 0) then
             do k = 1, nlayer
                layerSnowDensity_final(ij,k) = 0.
                snowIce_final(ij,k) = 0.
                snowLiq_final(ij,k) = 0.
                snowTProf_final(ij,k) = tm ! Default in JULES
                layerSnowGrain_final(ij,k) = 50. ! Default in JULES
                layerSnowDepth_final(ij,k) = 0.
             end do
             actSnowNL_final(ij) = 0
             grndSnow_final(ij) = 0. ! Always zero in PS41
             SWE_final(ij) = 0.
             snowDensity_final(ij) = 109. ! Default in JULES
             snowDepth_final(ij) = 0.
             snowGrain_final(ij) = 50. ! Default in JULES
             surftSnow_final(ij) = 0.
          else

             ! We have snow, so finish relayering and copy to final arrays
             call step_10()
             call step_11()
             do k = 1, nlayer
                layerSnowDensity_final(ij,k) = new_rho_snow(k)
                snowIce_final(ij,k) = new_sice(k)
                snowLiq_final(ij,k) = new_sliq(k)
                snowTProf_final(ij,k) = new_tsnow(k)
                layerSnowGrain_final(ij,k) = new_rgrainl(k)
                layerSnowDepth_final(ij,k) = new_ds(k)
             end do ! k
             actSnowNL_final(ij) = new_nsnow
             grndSnow_final(ij) = 0. ! Always zero in PS41
             SWE_final(ij) = new_snowmass
             snowDensity_final(ij) = new_rho_grnd
             snowDepth_final(ij) = new_snowdepth
             snowGrain_final(ij) = new_rgrain
             surftSnow_final(ij) = new_snowmass

          end if

       end do ! c
    end do ! r

    ! Free up internal memory, but keep "final" arrays for subsequent
    ! copying to external routine.
    call deallocate_step_arrays()
    call nullify_pointers()

  end subroutine LVT_proc_jules_ps41_ens_snow

  ! Fetch appropriate "final" array based on requested variable name.
  ! We also deallocate single-layer "final" array after copying to
  ! reduce memory usage; the multi-layer arrays must remain allocated,
  ! since only a single vertical slice is copied back by this routine.
  subroutine LVT_fetch_jules_ps41_ens_snow_final(nc, nr, data, k, short_name, &
       is_ps41_snow_var)

    ! Defaults
    implicit none

    ! Arguments
    integer, intent(in) :: nc
    integer, intent(in) :: nr
    real, intent(inout) :: data(nc*nr)
    integer, intent(in) :: k
    character(len=*), intent(in) :: short_name
    logical, intent(out) :: is_ps41_snow_var

    is_ps41_snow_var = .true. ! First guess

    select case(trim(short_name))
    case ("ActSnowNL_inst")
       data(:) = actSnowNL_final(:)
       deallocate(actSnowNL_final)
    case ("GrndSnow_inst")
       data(:) = grndSnow_final(:)
       deallocate(grndSnow_final)
    case ("LayerSnowDensity_inst")
       data(:) = layerSnowDensity_final(:,k)
    case ("LayerSnowDepth_inst")
       data(:) = layerSnowDepth_final(:,k)
    case ("LayerSnowGrain_inst")
       data(:) = layerSnowGrain_final(:,k)
    case ("SWE_inst")
       data(:) = SWE_final(:)
       deallocate(SWE_final)
    case ("SnowDensity_inst")
       data(:) = snowDensity_final(:)
       deallocate(snowDensity_final)
    case ("SnowDepth_inst")
       data(:) = snowDepth_final(:)
       deallocate(snowDepth_final)
    case ("SnowGrain_inst")
       data(:) = snowGrain_final(:)
       deallocate(snowGrain_final)
    case ("SnowIce_inst")
       data(:) = snowIce_final(:,k)
    case ("SnowLiq_inst")
       data(:) = snowLiq_final(:,k)
    case ("SnowTProf_inst")
       data(:) = SnowTProf_final(:,k)
    case ("SurftSnow_inst")
       data(:) = surftSnow_final(:)
       deallocate(surftSnow_final)
    case default
       is_ps41_snow_var = .false.
    end select

  end subroutine LVT_fetch_jules_ps41_ens_snow_final

  ! Deallocate arrays for final data
  subroutine LVT_cleanup_jules_ps41_ens_snow()
    implicit none
    if (allocated(actSnowNL_final)) deallocate(actSnowNL_final)
    if (allocated(grndSnow_final)) deallocate(grndSnow_final)
    if (allocated(layerSnowDensity_final)) deallocate(layerSnowDensity_final)
    if (allocated(layerSnowDepth_final)) deallocate(layerSnowDepth_final)
    if (allocated(layerSnowGrain_final)) deallocate(layerSnowGrain_final)
    if (allocated(SWE_final)) deallocate(SWE_final)
    if (allocated(snowDensity_final)) deallocate(snowDensity_final)
    if (allocated(snowDepth_final)) deallocate(snowDepth_final)
    if (allocated(snowGrain_final)) deallocate(snowGrain_final)
    if (allocated(snowIce_final)) deallocate(snowIce_final)
    if (allocated(snowLiq_final)) deallocate(snowLiq_final)
    if (allocated(snowTProf_final)) deallocate(snowTProf_final)
    if (allocated(surftSnow_final)) deallocate(surftSnow_final)
  end subroutine LVT_cleanup_jules_ps41_ens_snow

  ! Set metadata structure for SWE
  subroutine LVT_set_SWE_metadata(SWE)
    use LVT_LISoutputHandlerMod, only: LVT_lismetadataEntry
    implicit none
    type(LVT_lismetadataEntry), intent(inout) :: SWE
    SWE%short_name = "SWE"
    SWE%long_name = "snow water equivalent"
    SWE%standard_name = "liquid_water_content_of_surface_snow"
    SWE%units = "kg m-2"
    SWE%nunits = 1
    SWE%format = 'F'
    SWE%vlevels = 1
    SWE%timeAvgOpt = 0
    allocate(SWE%unittypes(1))
    SWE%unittypes(1) = "kg m-2"
    allocate(SWE%valid_min(1))
    SWE%valid_min = (/ 0.0 /)
    allocate(SWE%valid_max(1))
    SWE%valid_max = (/ 2000.0 /)
    SWE%varid_def = -99
    SWE%selectOpt = 1
  end subroutine LVT_set_SWE_metadata

  ! Set metadata structure for SnowDensity
  subroutine LVT_set_SnowDensity_metadata(SnowDensity)
    use LVT_LISoutputHandlerMod, only: LVT_lismetadataEntry
    implicit none
    type(LVT_lismetadataEntry), intent(inout) :: SnowDensity
    SnowDensity%short_name = "SnowDensity"
    SnowDensity%long_name = "snowpack bulk density"
    SnowDensity%standard_name = "snowpack_bulk_density"
    SnowDensity%units = "kg m-3"
    SnowDensity%nunits = 1
    SnowDensity%format = 'F'
    SnowDensity%vlevels = 1
    SnowDensity%timeAvgOpt = 0
    allocate(SnowDensity%unittypes(1))
    SnowDensity%unittypes(1) = "kg m-3"
    allocate(SnowDensity%valid_min(1))
    SnowDensity%valid_min = (/ 109.0 /)
    allocate( SnowDensity%valid_max(1))
    SnowDensity%valid_max = (/ 917.0 /)
    SnowDensity%varid_def = -99
    SnowDensity%selectOpt = 1
  end subroutine LVT_set_SnowDensity_metadata

  ! Set metadata structure for LayerSnowDensity
  subroutine LVT_set_LayerSnowDensity_metadata(LayerSnowDensity)
    use LVT_LISoutputHandlerMod, only: LVT_lismetadataEntry
    implicit none
    type(LVT_lismetadataEntry), intent(inout) :: LayerSnowDensity
    LayerSnowDensity%short_name = "LayerSnowDensity"
    LayerSnowDensity%long_name = "snow_density_for_each_layer"
    LayerSnowDensity%standard_name = "snow_density_for_each_layer"
    LayerSnowDensity%units = "kg m-3"
    LayerSnowDensity%nunits = 1
    LayerSnowDensity%format = 'F'
    LayerSnowDensity%vlevels = 3
    LayerSnowDensity%timeAvgOpt = 0
    allocate(LayerSnowDensity%unittypes(1))
    LayerSnowDensity%unittypes(1) = "kg m-3"
    allocate(LayerSnowDensity%valid_min(1))
    LayerSnowDensity%valid_min = (/ 109.0 /)
    allocate(LayerSnowDensity%valid_max(1))
    LayerSnowDensity%valid_max = (/ 917.0 /)
    LayerSnowDensity%varid_def = -99
    LayerSnowDensity%selectOpt = 1
  end subroutine LVT_set_LayerSnowDensity_metadata

  ! Set metadata structure for SnowGrain
  subroutine LVT_set_SnowGrain_metadata(SnowGrain)
    use LVT_LISoutputHandlerMod, only: LVT_lismetadataEntry
    implicit none
    type(LVT_lismetadataEntry), intent(inout) :: SnowGrain
    SnowGrain%short_name = "SnowGrain"
    SnowGrain%long_name = "snow grain size"
    SnowGrain%standard_name = "snow_grain_size"
    SnowGrain%units = "micron"
    SnowGrain%nunits = 1
    SnowGrain%format = 'F'
    SnowGrain%vlevels = 1
    SnowGrain%timeAvgOpt = 0
    allocate(SnowGrain%unittypes(1))
    SnowGrain%unittypes(1) = "micron"
    allocate(SnowGrain%valid_min(1))
    SnowGrain%valid_min = (/ 50.0 /)
    allocate(SnowGrain%valid_max(1))
    SnowGrain%valid_max = (/ 2000.0 /)
    SnowGrain%varid_def = -99
    SnowGrain%selectOpt = 1
  end subroutine LVT_set_SnowGrain_metadata

  ! Set metadata structure for SnowDepth
  subroutine LVT_set_SnowDepth_metadata(SnowDepth)
    use LVT_LISoutputHandlerMod, only: LVT_lismetadataEntry
    implicit none
    type(LVT_lismetadataEntry), intent(inout) :: SnowDepth
    SnowDepth%short_name = "SnowDepth"
    SnowDepth%long_name = "snow depth"
    SnowDepth%standard_name = "snow_depth"
    SnowDepth%units = "m"
    SnowDepth%nunits = 1
    SnowDepth%format = 'F'
    SnowDepth%vlevels = 1
    SnowDepth%timeAvgOpt = 0
    allocate(SnowDepth%unittypes(1))
    SnowDepth%unittypes(1) = "m"
    allocate(SnowDepth%valid_min(1))
    SnowDepth%valid_min = (/ 0.0 /)
    allocate(SnowDepth%valid_max(1))
    SnowDepth%valid_max = (/ 100.0 /)
    SnowDepth%varid_def = -99
    SnowDepth%selectOpt = 1
  end subroutine LVT_set_SnowDepth_metadata

  ! Set metadata structure for GrndSnow
  subroutine LVT_set_grndsnow_metadata(grndSnow)
    use LVT_LISoutputHandlerMod, only: LVT_lismetadataEntry
    implicit none
    type(LVT_lismetadataEntry), intent(inout) :: grndSnow
    grndSnow%short_name = "GrndSnow"
    grndSnow%long_name = "snow on ground (beneath canopy)"
    grndSnow%standard_name = "snow_on_ground_beneath_canopy"
    grndSnow%units = "kg m-2"
    grndSnow%nunits = 1
    grndSnow%format = 'F'
    grndSnow%vlevels = 1
    grndSnow%timeAvgOpt = 0
    allocate(grndSnow%unittypes(1))
    grndSnow%unittypes(1) = "kg m-2"
    allocate(grndSnow%valid_min(1))
    grndSnow%valid_min = (/ 0.0 /)
    allocate(grndSnow%valid_max(1))
    grndSnow%valid_max = (/ 1000.0 /)
    grndSnow%varid_def = -99
    grndSnow%selectOpt = 1
  end subroutine LVT_set_grndsnow_metadata

  ! Set metadata structure for SurftSnow
  subroutine LVT_set_surftsnow_metadata(surftSnow)
    use LVT_LISoutputHandlerMod, only: LVT_lismetadataEntry
    implicit none
    type(LVT_lismetadataEntry), intent(inout) :: surftSnow
    surftSnow%short_name = "SurftSnow"
    surftSnow%long_name = "snow amount on tile"
    surftSnow%standard_name = "snow_amount_on_tile"
    surftSnow%units = "kg m-2"
    surftSnow%nunits = 1
    surftSnow%format = 'F'
    surftSnow%vlevels = 1
    surftSnow%timeAvgOpt = 0
    allocate(surftSnow%unittypes(1))
    surftSnow%unittypes(1) = "kg m-2"
    allocate(surftSnow%valid_min(1))
    surftSnow%valid_min = (/ 0.0 /)
    allocate(surftSnow%valid_max(1))
    surftSnow%valid_max = (/ 10000.0 /)
    surftSnow%varid_def = -99
    surftSnow%selectOpt = 1
  end subroutine LVT_set_surftsnow_metadata

  ! Deallocates internal arrays within structure
  subroutine LVT_deallocate_metadata(var)
    use LVT_LISoutputHandlerMod, only: LVT_lismetadataEntry
    implicit none
    type(LVT_lismetadataEntry), intent(inout) :: var
    if (allocated(var%unittypes)) deallocate(var%unittypes)
    if (allocated(var%valid_min)) deallocate(var%valid_min)
    if (allocated(var%valid_max)) deallocate(var%valid_max)
  end subroutine LVT_deallocate_metadata

end module LVT_557post_ps41_snowMod

