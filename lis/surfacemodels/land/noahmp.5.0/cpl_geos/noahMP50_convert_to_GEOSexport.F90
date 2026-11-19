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
  associate( nm => NoahMP50_struc(n)%noahmp50(t) )

        fveg = nm%fveg

        ! --- Total Runoff flux ---
        nm%totalrunoff = (nm%runsf + nm%runsb) / nm%soiltimestep

        ! --- Heat and moisture exchange coefficients ---
        nm%ch_f = nm%rhoair * nm%ch

        ! --- Momentum exchange coefficient ---
        nm%cm_f = nm%rhoair * nm%cm

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
           nm%soilmoist_rz = rznumer / rzdenom
        else
           nm%soilmoist_rz = LIS_rc%udef
        endif

        ! --- canopy interception water storage ---
        nm%canopint = nm%canliq + nm%canice

        ! --- snow depth, snow-covered-area form ---
        snc = nm%snowc 
        if (snc > EPS) then
           nm%snowdepth_sca = nm%snowh / snc
        else
           nm%snowdepth_sca = nm%snowh
        endif

        ! --- band-resolved albedo ---
        nm%albedo_visdir = nm%albsfcdir(1)   
        nm%albedo_nirdir = nm%albsfcdir(2)   
        nm%albedo_visdif = nm%albsfcdif(1)   
        nm%albedo_nirdif = nm%albsfcdif(2)   

        !---- Vegetation area index ---
        nm%vai = nm%LAI + nm%SAI

        !---- Friction velocity ---
        nm%ustar = fveg*nm%ustarveg + (1.0-fveg)*nm%ustarbare        
    end associate
  enddo

end subroutine NoahMP50_convert_to_GEOSexport
