module NoahMP50_geosCplMod
! Storage for GEOS-coupling CONVERSION/DERIVED variables.
! These are NOT produced by Noah-MP physics as such; they are synthesised
! only on the GEOS-coupled path. Standalone runs never allocate these.
  use LIS_coreMod, only : LIS_rc
  implicit none
  private
  public :: noahmp50_geoscpl_dec, NoahMP50_geoscpl
  public :: NoahMP50_geoscpl_init, NoahMP50_geoscpl_finalize

  integer, parameter :: NBAND = 2   ! 1 = visible, 2 = near-infrared

  type :: noahmp50_geoscpl_dec
     real          :: totalrunoff    ! sfcrunoff + udrrunoff        [mm/soil_ts]
     real          :: ch_f           ! heat exchange coeff          [kg m-2 s-1]
     real          :: cm_f           ! gridcell momentum exchange coeff [kg m-2 s-1] 
     real          :: soilmoist_rz   ! depth-weighted root-zone SM  [m3/m3]
     real          :: canopint       ! canliq + canice              [kg/m2]
     real          :: snowdepth_sca  ! snowh/snowc (SCA form)       [m]
     real          :: roughness_f    ! z0 (also feeds Z0H)          [m]
     real          :: lai            ! LAI + SAI                    [m2/m2]
     real          :: ustar          ! fveg-weighted friction velocity [m/s]
     ! band-resolved albedo (1=vis, 2=nir); harvested from NoahmpIO in main
     real, pointer :: albedo_dir(:)  ! (NBAND): AlbedoVisDir, AlbedoNIRDir
     real, pointer :: albedo_dif(:)  ! (NBAND): AlbedoVisDif, AlbedoNIRDif
  end type noahmp50_geoscpl_dec

  type :: noahmp50_geoscpl_type
     type(noahmp50_geoscpl_dec), allocatable :: tile(:)
  end type noahmp50_geoscpl_type

  type(noahmp50_geoscpl_type), allocatable :: NoahMP50_geoscpl(:)

contains

  subroutine NoahMP50_geoscpl_init(n)
    ! Idempotent: safe to call every coupling step. Only allocates this
    ! nest's tile(:) store the first time it is invoked.
    integer, intent(in) :: n
    integer :: t, npatch
    if (.not. allocated(NoahMP50_geoscpl)) &
         allocate(NoahMP50_geoscpl(LIS_rc%nnest))
    if (allocated(NoahMP50_geoscpl(n)%tile)) return
    npatch = LIS_rc%npatch(n, LIS_rc%lsm_index)
    allocate(NoahMP50_geoscpl(n)%tile(npatch))
    do t = 1, npatch
       allocate(NoahMP50_geoscpl(n)%tile(t)%albedo_dir(NBAND))
       allocate(NoahMP50_geoscpl(n)%tile(t)%albedo_dif(NBAND))
       NoahMP50_geoscpl(n)%tile(t)%albedo_dir = LIS_rc%udef
       NoahMP50_geoscpl(n)%tile(t)%albedo_dif = LIS_rc%udef
    enddo
  end subroutine NoahMP50_geoscpl_init

  subroutine NoahMP50_geoscpl_finalize(n)
    integer, intent(in) :: n
    integer :: t
    if (.not. allocated(NoahMP50_geoscpl)) return
    do t = 1, size(NoahMP50_geoscpl(n)%tile)
       if (associated(NoahMP50_geoscpl(n)%tile(t)%albedo_dir)) &
            deallocate(NoahMP50_geoscpl(n)%tile(t)%albedo_dir)
       if (associated(NoahMP50_geoscpl(n)%tile(t)%albedo_dif)) &
            deallocate(NoahMP50_geoscpl(n)%tile(t)%albedo_dif)
    enddo
    deallocate(NoahMP50_geoscpl(n)%tile)
  end subroutine NoahMP50_geoscpl_finalize

end module NoahMP50_geosCplMod
