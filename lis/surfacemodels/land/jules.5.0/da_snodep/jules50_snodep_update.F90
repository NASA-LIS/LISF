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
! !ROUTINE: jules50_snodep_update
! \label{jules50_snodep_update}
!
! !REVISION HISTORY:
!  18 Jun 2019: Yeosang Yoon; Initial specification
!
! !INTERFACE

subroutine jules50_snodep_update(n, t, dsneqv, dsnowh)

  use LIS_coreMod
  use jules50_lsmMod
  USE jules_snow_mod, ONLY:                                                     &
    nsmax,                                                                      &
      ! Maximum possible number of snow layers.
    r0                                                                         
      ! Grain size for fresh snow (microns).

  USE compactsnow_mod,     ONLY: compactsnow  
  USE relayersnow_mod,     ONLY: relayersnow
  USE snowgrain_mod,       ONLY: snowgrain
  USE layersnow_mod,       ONLY: layersnow
  USE ancil_info,          ONLY: land_pts, nsurft, nsoilt
  USE c_0_dg_c,            ONLY: tm
  USE jules_radiation_mod, ONLY: l_snow_albedo, l_embedded_snow


  implicit none
!
! !DESCRIPTION:
!  This subroutine updates relevant snow prognostics based
!  on the update to the total SWE (dsneqv) and total
!  snow depth (dsnowh). The updated variables include
!  number of snow layers, snice, snliq, snow temperature
!  snow grain size, and snow thickness.
!
! !ARGUMENTS:
  integer, intent(in)  :: n
  integer, intent(in)  :: t
  real                 :: dsneqv !mm
  real                 :: dsnowh !m
!EOP

  integer         :: p, q, i, j, k, n_can
  integer         :: surft_pts(nsurft)            ! Number of tile points.
  integer         :: surft_index(land_pts,nsurft) ! Index of tile points.
  integer         :: nsnow(land_pts,nsurft)       ! Number of snow layers.
  real            :: timestep                     ! Timestep length (s).
  

  real            :: sice(land_pts,nsurft,nsmax)     ! Ice content of snow layers (kg/m2).
  real            :: sliq(land_pts,nsurft,nsmax)     ! Liquid content of snow layers (kg/m2).
  real            :: rgrain(land_pts,nsurft)         ! Snow surface grain size (microns).
  real            :: rgrainl(land_pts,nsurft,nsmax)  ! Snow layer grain size (microns).
  real            :: rho_snow_grnd(land_pts,nsurft)  ! Snowpack bulk density (kg/m3).
  real            :: snow_grnd(land_pts,nsurft)      ! Snow beneath canopy (kg/m2).
  real            :: snow_surft(land_pts,nsurft)     ! Snow mass (kg/m2).
  real            :: snowdepth(land_pts,nsurft)      ! Snow depth (m).
  real            :: tsnow(land_pts,nsurft,nsmax)    ! Snow layer temperatures (K).
  real            :: ds(land_pts,nsurft,nsmax)       ! Snow layer thicknesses (m).
  real            :: rho_snow(land_pts,nsurft,nsmax) ! Snow layer densities (kg/m3).  
  real            :: t_soil1_soilt(land_pts,nsoilt)  ! Soil surface layer temperature (K).
  real            :: tstar_surft(land_pts,nsurft)    ! Tile surface temperature (K).

!-----------------------------------------------------------------------------
! Local arrays
!-----------------------------------------------------------------------------
REAL ::                                                                       &
  rho0(land_pts),                                                             &
    ! Density of fresh snow (kg/m3).
    ! Where nsnow=0, rho0 is the density of the snowpack.
  snowfall(land_pts),                                                         &
    ! Total frozen precip reaching the ground in timestep
    ! (kg/m2) - includes any canopy unloading.
  snowmass(land_pts),                                                         &
    ! Snow mass on the ground (kg/m2).
  rgrain0(land_pts),                                                          &
    ! Fresh snow grain size (microns).
  sice0(land_pts),                                                            &
   ! Ice content of fresh snow (kg/m2).
   ! Where nsnow=0, sice0 is the mass of the snowpack.
  tsnow0(land_pts)
    ! Temperature of fresh snow (K).
!-----------------------------------------------------------------------------

  !---------------------------------------------------------------------------
  ! Get variables
  !---------------------------------------------------------------------------
  timestep                       = jules50_struc(n)%ts
  surft_pts(nsurft)              = jules50_struc(n)%jules50(t)%tile_pts(nsurft)
  surft_index(land_pts,nsurft)   = jules50_struc(n)%jules50(t)%tile_index(nsurft)
  nsnow(land_pts,nsurft)         = jules50_struc(n)%jules50(t)%nsnow(nsurft)
  snow_surft(land_pts,nsurft)    = jules50_struc(n)%jules50(t)%snow_tile(nsurft)   ! ntiles, Lying snow on tiles (kg m-2)
  snow_grnd(land_pts,nsurft)     = jules50_struc(n)%jules50(t)%snow_grnd(nsurft)   ! ntiles, Snow on the ground (kg m-2)
  t_soil1_soilt(land_pts,nsoilt) = jules50_struc(n)%jules50(t)%t_soil(nsoilt)      ! sm_levels, Sub-surface temperatures (K)
  rho_snow_grnd(land_pts,nsurft) = jules50_struc(n)%jules50(t)%rho_snow_grnd(nsurft)  
  tstar_surft(land_pts,nsurft)   = jules50_struc(n)%jules50(t)%tstar_tile(nsurft)  ! ntiles, Tile surface temperatures (K)
  rgrain(land_pts,nsurft)        = jules50_struc(n)%jules50(t)%rgrain(nsurft)      ! snow surface grain size
  snowdepth(land_pts,nsurft)     = jules50_struc(n)%jules50(t)%snowdepth(nsurft)
  snowmass(land_pts)             = jules50_struc(n)%jules50(t)%snow_mass_ij
   
  DO i=1, nsmax
    sice(land_pts,nsurft,i)     = jules50_struc(n)%jules50(t)%sice(nsurft,i)      ! snow ice
    sliq(land_pts,nsurft,i)     = jules50_struc(n)%jules50(t)%sliq(nsurft,i)      ! snow liquid water
    tsnow(land_pts,nsurft,i)    = jules50_struc(n)%jules50(t)%tsnow(nsurft,i)     ! snow temp.
    rho_snow(land_pts,nsurft,i) = jules50_struc(n)%jules50(t)%rho_snow(nsurft,i)  
    rgrainl(land_pts,nsurft,i)  = jules50_struc(n)%jules50(t)%rgrainl(nsurft,i)   ! snow grain size
    ds(land_pts,nsurft,i)       = jules50_struc(n)%jules50(t)%ds(nsurft,i)        ! snow layer thicknesses
  END DO

DO p = 1,nsurft

  DO k = 1,surft_pts(p)
    i = surft_index(k,p)
    ! Update snow depth and swe from observations
    snowdepth(i,p) = snowdepth(i,p) + dsnowh
    snowmass(i) = snowmass(i) + dsneqv

    IF ( snowdepth(i,p) < 1.0E-6 .OR. snowmass(i) < 1.0E-3 ) THEN
      snowdepth(i,p)  = 0.0
      snowmass(i)     = 0.0
      snowfall        = 0.0
      nsnow(i,p)      = 0
      ds(:,p,:)       = 0.0
      sice(:,p,:)     = 0.0
      sliq(:,p,:)     = 0.0
      rgrain(i,p)     = r0
      rgrainl(:,p,:)  = r0
      rho_snow(:,p,:) = 0.0
      tsnow(:,p,:)    = tm
    END IF
 
    IF ( snowdepth(i,p) > 0 ) THEN 
      rho_snow_grnd(i,p) = snowmass(i) / snowdepth(i,p)
    END IF
  END DO

  !---------------------------------------------------------------------------
  ! Divide snow pack into layers
  !---------------------------------------------------------------------------
  CALL layersnow ( land_pts, surft_pts(p), surft_index(:, p),                 &
                   snowdepth(:,p), nsnow(:,p), ds(:,p,:) )

  snowfall = 0.0
  !---------------------------------------------------------------------------
  ! Growth of snow grains
  !---------------------------------------------------------------------------
  IF ( l_snow_albedo .OR. l_embedded_snow ) THEN
    CALL snowgrain ( land_pts, surft_pts(p), timestep, nsnow(:,p),            &
                     surft_index(:,p), sice(:,p,:), snowfall,                 &
                     snowmass, tsnow(:,p,:), tstar_surft(:,p),                &
                     rgrain(:,p), rgrainl(:,p,:), rgrain0 )
  ELSE
    ! Default initialization required for bit-comparison in the UM.
    rgrain0(:) = r0
  END IF

  IF ( nsmax > 0 ) THEN
    !---------------------------------------------------------------------------
    ! Snow thermodynamics and hydrology
    !---------------------------------------------------------------------------
    DO k = 1,surft_pts(p)
      i = surft_index(k,p)

      IF ( nsnow(i,p) == 0 ) THEN        
        ! Diagnose snow depth.
        snowdepth(i,p) = snowmass(i) / rho_snow_grnd(i,p)

        ! Set values for the surface layer.
        rho0(i)   = rho_snow_grnd(i,p)
        sice0(i)  = snowmass(i)
        tsnow0(i) = MIN( t_soil1_soilt(i,p), tm )
      ELSE
        ! Diagnose layer densities
        DO q = 1,nsnow(i,p)
          IF ( ds(i,p,q) > EPSILON(ds) ) THEN
            rho_snow(i,p,q) = (sice(i,p,q) + sliq(i,p,q)) / ds(i,p,q)
          END IF
        END DO

        rho0(i)   = rho_snow_grnd(i,p)
        sice0(i)  = dsneqv
        tsnow0(i) = tsnow(i,p,1)

        ! Diagnose total snow depth and mass
        snowdepth(i,p) = sice0(i) / rho0(i)
        snowmass(i)  = sice0(i)
        DO q = 1,nsnow(i,p)
          snowdepth(i,p) = snowdepth(i,p) + ds(i,p,q)
          snowmass(i)    = snowmass(i) + sice(i,p,q) + sliq(i,p,q)
        END DO
      END IF
    END DO

    !-------------------------------------------------------------------------
    ! Mechanical compaction of snow
    !-------------------------------------------------------------------------
    CALL compactsnow ( land_pts, surft_pts(p), timestep, nsnow(:,p),          &
                       surft_index(:,p), sice(:,p,:), sliq(:,p,:),            &
                       tsnow(:,p,:), rho_snow(:,p,:), ds(:,p,:) )

    !-------------------------------------------------------------------------
    ! Redivide snowpack after changes in depth, conserving mass and energy
    !-------------------------------------------------------------------------
    CALL relayersnow ( land_pts, surft_pts(p), surft_index(:,p),              &
                       rgrain0, rho0, sice0, snowfall,                        &
                       snowmass, tsnow0, nsnow(:,p), ds(:,p,:),               &
                       rgrain(:,p), rgrainl(:,p,:), sice(:,p,:),              &
                       rho_snow_grnd(:,p), sliq(:,p,:),                       &
                       tsnow(:,p,:), rho_snow(:,p,:),                         &
                       snowdepth(:,p) )

  END IF  !  NSMAX>0

  !---------------------------------------------------------------------------
  ! Copy into final snow mass variables
  !---------------------------------------------------------------------------
  DO k = 1,surft_pts(p)
    i = surft_index(k,p)
    snow_surft(i,p) = snowmass(i)
  END DO

  !---------------------------------------------------------------------------
  ! Update snow variables
  !---------------------------------------------------------------------------
  DO k = 1,surft_pts(p)
    i = surft_index(k,p)
    IF (nsnow(i,p) > 0) THEN
      snowdepth(i,p) = 0.0
      snowmass(i)    = 0.0
      DO q = 1,nsnow(i,p)
        snowdepth(i,p) = snowdepth(i,p) + ds(i,p,q)
        snowmass(i)    = snowmass(i) + sice(i,p,q) + sliq(i,p,q)
      END DO
    END IF

    IF ( snowdepth(i,p) < 1.0E-6 .OR. snowmass(i) < 1.0E-3 ) THEN
      snowdepth(i,p)  = 0.0
      snowmass(i)     = 0.0
      snowfall        = 0.0
      nsnow(i,p)      = 0
      ds(:,p,:)       = 0.0
      sice(:,p,:)     = 0.0
      sliq(:,p,:)     = 0.0
      rgrain(i,p)     = r0
      rgrainl(:,p,:)  = r0
      rho_snow(:,p,:) = 0.0
      tsnow(:,p,:)    = tm
    END IF

    jules50_struc(n)%jules50(t)%snowdepth(p)     = snowdepth(i,p)
    jules50_struc(n)%jules50(t)%snow_mass_ij     = snowmass(i)
  END DO

  DO k = 1,surft_pts(p)
    i = surft_index(k,p)
    DO q=1, nsmax
      jules50_struc(n)%jules50(t)%sice(p,q)        = sice(i,p,q)
      jules50_struc(n)%jules50(t)%sliq(p,q)        = sliq(i,p,q)
      jules50_struc(n)%jules50(t)%tsnow(p,q)       = tsnow(i,p,q)
      jules50_struc(n)%jules50(t)%rgrainl(p,q)     = rgrainl(i,p,q)
      jules50_struc(n)%jules50(t)%ds(p,q)          = ds(i,p,q)
      jules50_struc(n)%jules50(t)%rho_snow(p,q)    = rho_snow(i,p,q)
    END DO

    jules50_struc(n)%jules50(t)%rgrain(p)        = rgrain(i,p)
    jules50_struc(n)%jules50(t)%rho_snow_grnd(p) = rho_snow_grnd(i,p)
    jules50_struc(n)%jules50(t)%snow_tile(p)     = snow_surft(i,p)
    jules50_struc(n)%jules50(t)%nsnow(p)         = nsnow(i,p)
  END DO

END DO  !  tiles

end subroutine jules50_snodep_update

