!BOP
! !ROUTINE: NoahMP50_convert_to_GEOSexport
! !DESCRIPTION:
!  Converts native Noah-MP state into GEOS/CLSM-ready export variables.
!  Invoked only on the GEOS-coupled path, after NoahMP50_main.
!EOP
subroutine NoahMP50_convert_to_GEOSexport(n)
! !USES:
  use LIS_coreMod
  use LIS_logMod
  use NoahMP50_lsmMod
  use NoahMP50_geosCplMod
  use NoahmpIOVarType        ! for NoahmpIO%dzs(:)

  implicit none
! !ARGUMENTS:
  integer, intent(in) :: n
! !LOCAL:
  integer :: t, i, nrz
  real    :: fveg, snc, cumdepth, rzdenom, rznumer
  real, parameter :: EPS     = 1.0e-6
  real, parameter :: RZ_TGT  = 1.0     ! root-zone target depth [m]

  do t = 1, LIS_rc%npatch(n, LIS_rc%lsm_index)
     associate( nm => NoahMP50_struc(n)%noahmp50(t), &
                gx => NoahMP50_geoscpl(n)%tile(t) )

        fveg = nm%fveg

        ! --- total runoff ---
        ! This is total runoff, but I think geos needs runoff flux 
        ! gx%totalrunoff = nm%sfcrunoff + nm%udrrunoff 
        
        ! --- Total Runoff flux [unit = kg m-2 s-1] ---
        ! NoahMP runsf and runsb are timestep runoff amounts [mm]
        ! ts is the land model time step
        gx%totalrunoff = (nm%runsf + nm%runsb) / NoahMP50_struc(n)%soiltimestep

        ! --- Heat and moisture exchange coefficients  ---
        ! NoahMP ch is in m/s
        ! GEOS expects kg m-2 s-1
        ! air density (rhoair in kg/m3) 
        ! This is heat coefficient. Assume heat and moisture coefficients are equal
           gx%ch_f = nm%rhoair * nm%ch

        ! --- Momentum exchange coefficient (unit = kg m-2 s-1) ---
        ! Noahmp cm is gridmean in m/s
        ! Multipy by air density for kg m-2 s-1 
          gx%cm_f = rm%rhoair * nm%cm

        ! --- root-zone soil moisture: depth-weighted mean over the top
        !     layers whose cumulative thickness first reaches/exceeds 1 m.
        !     Reads thicknesses from NoahmpIO%dzs so it tracks any retuning. ---
        !     unit = m3/m3
        cumdepth = 0.0 ; rznumer = 0.0 ; rzdenom = 0.0 ; nrz = 0
        do i = 1, size(NoahmpIO%dzs)
           rznumer  = rznumer + NoahmpIO%dzs(i)*nm%smc(i)
           rzdenom  = rzdenom + NoahmpIO%dzs(i)
           cumdepth = cumdepth + NoahmpIO%dzs(i)
           nrz = i
           if (cumdepth >= RZ_TGT - EPS) exit    ! include boundary-crossing layer
        enddo
        if (rzdenom > EPS) then
           gx%soilmoist_rz = rznumer / rzdenom
        else
           gx%soilmoist_rz = LIS_rc%udef
        endif


        ! --- canopy interception water storage (unit = kg/m2) ---
        ! NoahMP canliq and canice are in mm water equivalent
        ! 1 mm water equiv = 1 kg/m2, so no conversion needed for units 
        gx%canopint = nm%canliq + nm%canice

        ! --- snow depth, snow-covered-area form (guard snowc=0; unit=m) ---
        snc = nm%snowc 
        if (snc > EPS) then
           gx%snowdepth_sca = nm%snowh / snc
        else
           gx%snowdepth_sca = nm%snowh
        endif

        ! --- band-resolved albedo: direct SELECT from harvested native
        !     arrays (1=vis, 2=nir). Requires nm%albd(:)/nm%albi(:) to be
        !     populated in NoahMP50_main (see note). ---
        ! unit = dimensionless 
        gx%albedo_dir(1) = nm%albsfcdir(1)   ! AlbedoVisDir
        gx%albedo_dir(2) = nm%albsfcdir(2)   ! AlbedoNIRDir
        gx%albedo_dif(1) = nm%albsfcdif(1)   ! AlbedoVisDif
        gx%albedo_dif(2) = nm%albsfcdif(2)   ! AlbedoNIRDif

        !---- Vegetation area index (unit = m2/m2)--
        ! GEOS LAI = NoahMP leaf area index + stem area index
        gx%lai = nm%LAI + nm%SAI

        !---- Friction velocity (units = m/s)---
        ! GEOS USTAR = sum of fveg-weighted vegetated and bare ground
        ! USTAR = UstarVeg * VegFrac + UstarBare * (1 - VegFrac)
        gx%ustar = fveg*nm%ustarveg + (1.0-fveg)*nm%ustarbare        

     end associate
  enddo

end subroutine NoahMP50_convert_to_GEOSexport
