!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LIS_misc.h"
!BOP
!
! !ROUTINE: jules50_readrst
! \label{jules50_readrst}
!
! !REVISION HISTORY:
! 16 May 2016; Shugong Wang; initial implementation for JULES 5.0 
!
!
! !INTERFACE:
subroutine jules50_readrst()
! !USES:
    use LIS_coreMod, only    : LIS_rc, LIS_masterproc
    use LIS_historyMod, only : LIS_readvar_restart
    use LIS_logMod, only     : LIS_logunit, LIS_endrun, &
                               LIS_getNextUnitNumber,   &
                               LIS_releaseUnitNumber,   &
                               LIS_verify                
    use jules50_lsmMod

#if (defined USE_NETCDF3 || defined USE_NETCDF4)
    use netcdf
#endif

!
! !DESCRIPTION:
!  This program reads restart files for jules50.  This
!  includes all relevant water/energy storages and tile information.
!  The following is the list of variables specified in the jules50
!  restart file:
!  \begin{verbatim}
!    nc, nr, nsurfts             - grid and tile space dimensions
!    tsoil_deep                 - jules50 deep soil temperatures [K]
!    wood_prod_fast             - jules50 Fast-turnover wood product C pool [-]
!    wood_prod_med              - jules50 Medium-turnover wood product C pool [-]
!    wood_prod_slow             - jules50 Slow-turnover wood product C pool [-]
!    frac_agr_prev              - jules50 Agricultural fraction from previous TRIFFID call [-]
!    frac_agr                   - jules50 Agricultural fraction [-]
!    n_inorg                    - jules50 Gridbox Inorganic N pool [kg m-2]
!    canht_ft                   - jules50 Canopy height [m]
!    canopy                     - jules50 Surface/canopy water for snow-free land tiles [kg m-2]
!    ns                         - jules50 Soil Organic Nitrogen [kg m-2]
!    cs                         - jules50 Soil carbo [kg m-2]
!    gc                         - jules50 Stomatal conductance to evaporation for land tiles [m s-1]
!    gs                         - jules50 Stomatal conductance to evaporation [m s-1]
!    lai                        - jules50 LAI of plant functional types [-]
!    rgrain                     - jules50 Snow surface grain size on tiles [micron]
!    smc                        - jules50 Soil moisture in a layer at the surface [kg m-2]
!    smcl                       - jules50 Soil moisture content of layers [kg m-2]
!    snow_tile                  - jules50 Lying snow on tiles [kg m-2]
!    snow_grnd                  - jules50 Snow on the ground [kg m-2]
!    soot                       - jules50 Snow soot content (kg kg-1) [kg/kg]
!    t_soil                     - jules50 Sub-surface temperatures [K]
!    tstar_tile                 - jules50 Tile surface temperatures [K]
!    asteps_since_triffid       - jules50 Number of atmospheric timesteps since last call to TRIFFID [-]
!    lai_phen                   - jules50 LAI of PFTs after phenology. Required as separate variable for top-level argument list matching with VEG_IC2A [-]
!    c_veg                      - jules50 Total carbon content of the vegetation [kg m-2]
!    cv                         - jules50 Gridbox mean vegetation carbon [kg m-2]
!    fexp                       - jules50 Decay factor in Sat. Conductivity in water table layer [-]
!    ti_mean                    - jules50 Mean topographic index [-]
!    ti_sig                     - jules50 Standard deviation of topographic index [-]
!    zw                         - jules50 Water table depth [m]
!    sthzw                      - jules50 soil moist fraction in deep (water table) layer []
!    sthu                       - jules50 Unfrozen soil moisture content of the layers as a fraction of saturation [-]
!    sthf                       - jules50 Frozen soil moisture content of the layers as a fraction of saturation [-]
!    sice                       - jules50 Snow layer ice mass on tiles [kg m-2]
!    sliq                       - jules50 Snow layer liquid mass on tiles [kg m-2]
!    snowdepth                  - jules50 Snow depth on ground on tiles [m]
!    tsnow                      - jules50 Snow layer temperature [K]
!    rgrainl                    - jules50 Snow layer grain size on tiles [microns]
!    rho_snow_grnd              - jules50 Snowpack bulk density [kg/m3]
!    rho_snow                   - jules50 Snow layer densities [m]
!    ds                         - jules50 Snow layer thickness [m]
!  \end{verbatim}
!
!  The routines invoked are:
! \begin{description}
!   \item[LIS\_readvar\_restart](\ref{LIS_readvar_restart})\\
!      reads a variable from the restart file
!   \item[jules50\_coldstart](\ref{jules50_coldstart})\\
!      initializes the jules50 state variables
! \end{description}
!EOP
 
    implicit none
 
    integer           :: t, l, c, k 
    integer           :: nc, nr, npatch
    integer           :: n
    integer           :: ftn
    integer           :: status
    real, allocatable :: tmptilen(:)
    logical           :: file_exists
    integer           :: r_idx, c_idx 
    character*20      :: wformat
 
    do n=1, LIS_rc%nnest
        wformat = trim(JULES50_struc(n)%rformat)
        ! coldstart
        if(LIS_rc%startcode .eq. "coldstart") then  
            call jules50_coldstart(LIS_rc%lsm_index)
        ! restart
        elseif(LIS_rc%startcode .eq. "restart") then
            call jules50_initialize(LIS_rc%lsm_index)
            allocate(tmptilen(LIS_rc%npatch(n, LIS_rc%lsm_index)))
            ! check the existance of restart file
            inquire(file=JULES50_struc(n)%rfile, exist=file_exists)
            If (.not. file_exists) then 
                write(LIS_logunit,*) "jules50 restart file ", JULES50_struc(n)%rfile," does not exist "
                write(LIS_logunit,*) "Program stopping ..."
                call LIS_endrun
            endif
            write(LIS_logunit,*) "jules50 restart file used: ", JULES50_struc(n)%rfile
        
            ! open restart file
            if(wformat .eq. "binary") then
                ftn = LIS_getNextUnitNumber()
                open(ftn, file=JULES50_struc(n)%rfile, form="unformatted")
                read(ftn) nc, nr, npatch  !time, veg class, no. tiles
 
                ! check for grid space conflict
                if((nc .ne. LIS_rc%gnc(n)) .or. (nr .ne. LIS_rc%gnr(n))) then
                    write(LIS_logunit,*) JULES50_struc(n)%rfile, "grid space mismatch - jules50 halted"
                    call LIS_endrun
                endif
            
                if(npatch .ne. LIS_rc%glbnpatch_red(n, LIS_rc%lsm_index)) then
                    write(LIS_logunit,*) "restart tile space mismatch, halting..."
                    call LIS_endrun
                endif
            elseif(wformat .eq. "netcdf") then
#if (defined USE_NETCDF3 || defined USE_NETCDF4)
                status = nf90_open(path=JULES50_struc(n)%rfile, &
                                   mode=NF90_NOWRITE, ncid=ftn)
                call LIS_verify(status, "Error opening file "//JULES50_struc(n)%rfile)
#endif
            endif
 
            ! read: deep soil temperatures
            do l=1, JULES50_struc(n)%ns_deep ! TODO: check loop
                call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, tmptilen, varname="TSOIL_DEEP", &
                                         dim=l, vlevels = JULES50_struc(n)%ns_deep, wformat=wformat)
                do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
                    JULES50_struc(n)%jules50(t)%tsoil_deep(l) = tmptilen(t)
                enddo
            enddo
 
            ! read: Fast-turnover wood product C pool
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, JULES50_struc(n)%jules50%wood_prod_fast, &
                                     varname="WOOD_PROD_FAST", wformat=wformat)
 
            ! read: Medium-turnover wood product C pool
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, JULES50_struc(n)%jules50%wood_prod_med, &
                                     varname="WOOD_PROD_MED", wformat=wformat)
 
            ! read: Slow-turnover wood product C pool
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, JULES50_struc(n)%jules50%wood_prod_slow, &
                                     varname="WOOD_PROD_SLOW", wformat=wformat)
 
            ! read: Agricultural fraction from previous TRIFFID call
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, JULES50_struc(n)%jules50%frac_agr_prev, &
                                     varname="FRAC_AGR_PREV", wformat=wformat)
 
            ! read: Agricultural fraction
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, JULES50_struc(n)%jules50%frac_agr, &
                                     varname="FRAC_AGR", wformat=wformat)
 
            ! read: Gridbox Inorganic N pool
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, JULES50_struc(n)%jules50%n_inorg, &
                                     varname="N_INORG", wformat=wformat)
 
            ! read: Canopy height
            do l=1, JULES50_struc(n)%npft ! TODO: check loop
                call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, tmptilen, varname="CANHT_FT", &
                                         dim=l, vlevels = JULES50_struc(n)%npft, wformat=wformat)
                do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
                    JULES50_struc(n)%jules50(t)%canht_ft(l) = tmptilen(t)
                enddo
            enddo
 
            ! read: Surface/canopy water for snow-free land tiles
            do l=1, JULES50_struc(n)%nsurft ! TODO: check loop
                call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, tmptilen, varname="CANOPY", &
                                         dim=l, vlevels = JULES50_struc(n)%nsurft, wformat=wformat)
                do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
                    JULES50_struc(n)%jules50(t)%canopy(l) = tmptilen(t)
                enddo
            enddo
            
 
            ! read: Soil Organic Nitrogen
            do c=1, JULES50_struc(n)%dim_cslayer ! TODO: check loop
              do l=1, JULES50_struc(n)%dim_cs1 ! TODO: check loop
                  k = (c-1)*JULES50_struc(n)%dim_cslayer + l 
                  call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, tmptilen, varname="NS", &
                                           dim=k, vlevels =JULES50_struc(n)%dim_cslayer*JULES50_struc(n)%dim_cs1, wformat=wformat)
                  do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
                      !JULES50_struc(n)%jules50(t)%ns(l) = tmptilen(t)
                      JULES50_struc(n)%jules50(t)%ns(c, l) = tmptilen(t)
                  enddo
              enddo
            enddo
 
            ! read: Soil carbo
            do c=1, JULES50_struc(n)%dim_cslayer ! TODO: check loop
              do l=1, JULES50_struc(n)%dim_cs1 ! TODO: check loop
                  k = (c-1)*JULES50_struc(n)%dim_cslayer + l 
                  call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, tmptilen, varname="CS", &
                                           dim=k, vlevels =JULES50_struc(n)%dim_cslayer*JULES50_struc(n)%dim_cs1, wformat=wformat)
                  do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
                      JULES50_struc(n)%jules50(t)%cs_pool_soilt(c,l) = tmptilen(t)
                  enddo
              enddo
            enddo 
            ! read: Stomatal conductance to evaporation for land tiles
            do l=1, JULES50_struc(n)%nsurft ! TODO: check loop
                call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, tmptilen, varname="GC", &
                                         dim=l, vlevels = JULES50_struc(n)%nsurft, wformat=wformat)
                do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
                    JULES50_struc(n)%jules50(t)%gc_surft(l) = tmptilen(t)
                enddo
            enddo
 
            ! read: Stomatal conductance to evaporation
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, JULES50_struc(n)%jules50%gs_gb, &
                                     varname="GS", wformat=wformat)
 
            ! read: LAI of plant functional types
            do l=1, JULES50_struc(n)%npft ! TODO: check loop
                call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, tmptilen, varname="LAI", &
                                         dim=l, vlevels = JULES50_struc(n)%npft, wformat=wformat)
                do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
                    JULES50_struc(n)%jules50(t)%lai(l) = tmptilen(t)
                enddo
            enddo
 
            ! read: Snow surface grain size on tiles
            do l=1, JULES50_struc(n)%nsurft ! TODO: check loop
                call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, tmptilen, varname="RGRAIN", &
                                         dim=l, vlevels = JULES50_struc(n)%nsurft, wformat=wformat)
                do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
                    JULES50_struc(n)%jules50(t)%rgrain(l) = tmptilen(t)
                enddo
            enddo
 
            ! read: Soil moisture in a layer at the surface
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, JULES50_struc(n)%jules50%smc_soilt, &
                                     varname="SMC", wformat=wformat)
 
            ! read: Soil moisture content of layers
            do l=1, JULES50_struc(n)%sm_levels ! TODO: check loop
                call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, tmptilen, varname="SMCL", &
                                         dim=l, vlevels = JULES50_struc(n)%sm_levels, wformat=wformat)
                do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
                    JULES50_struc(n)%jules50(t)%smcl_soilt(l) = tmptilen(t)
                enddo
            enddo
 
            ! read: Lying snow on tiles
            do l=1, JULES50_struc(n)%nsurft ! TODO: check loop
                call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, tmptilen, varname="SNOW_TILE", &
                                         dim=l, vlevels = JULES50_struc(n)%nsurft, wformat=wformat)
                do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
                    JULES50_struc(n)%jules50(t)%snow_tile(l) = tmptilen(t)
                enddo
            enddo
 
            ! read: Snow on the ground
            do l=1, JULES50_struc(n)%nsurft ! TODO: check loop
                call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, tmptilen, varname="SNOW_GRND", &
                                         dim=l, vlevels = JULES50_struc(n)%nsurft, wformat=wformat)
                do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
                    JULES50_struc(n)%jules50(t)%snow_grnd(l) = tmptilen(t)
                enddo
            enddo
 
            ! read: Snow soot content (kg kg-1)
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, JULES50_struc(n)%jules50%soot_ij, &
                                     varname="SOOT", wformat=wformat)
 
            ! read: Sub-surface temperatures
            do l=1, JULES50_struc(n)%sm_levels ! TODO: check loop
                call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, tmptilen, varname="T_SOIL", &
                                         dim=l, vlevels = JULES50_struc(n)%sm_levels, wformat=wformat)
                do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
                    JULES50_struc(n)%jules50(t)%t_soil(l) = tmptilen(t)
                enddo
            enddo
 
            ! read: Tile surface temperatures
            do l=1, JULES50_struc(n)%nsurft ! TODO: check loop
                call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, tmptilen, varname="TSTAR_TILE", &
                                         dim=l, vlevels = JULES50_struc(n)%nsurft, wformat=wformat)
                do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
                    JULES50_struc(n)%jules50(t)%tstar_tile(l) = tmptilen(t)
                enddo
            enddo
 
            ! read: Number of atmospheric timesteps since last call to TRIFFID
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, JULES50_struc(n)%jules50%asteps_since_triffid, &
                                     varname="ASTEPS_SINCE_TRIFFID", wformat=wformat)
 
            ! read: LAI of PFTs after phenology. Required as separate variable for top-level argument list matching with VEG_IC2A
            do l=1, JULES50_struc(n)%npft ! TODO: check loop
                call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, tmptilen, varname="LAI_PHEN", &
                                         dim=l, vlevels = JULES50_struc(n)%npft, wformat=wformat)
                do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
                    JULES50_struc(n)%jules50(t)%lai_phen(l) = tmptilen(t)
                enddo
            enddo
 
            ! read: Total carbon content of the vegetation
            do l=1, JULES50_struc(n)%npft ! TODO: check loop
                call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, tmptilen, varname="C_VEG", &
                                         dim=l, vlevels = JULES50_struc(n)%npft, wformat=wformat)
                do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
                    JULES50_struc(n)%jules50(t)%c_veg(l) = tmptilen(t)
                enddo
            enddo
 
            ! read: Gridbox mean vegetation carbon
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, JULES50_struc(n)%jules50%cv, &
                                     varname="CV", wformat=wformat)
 
            ! read: Decay factor in Sat. Conductivity in water table layer
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, JULES50_struc(n)%jules50%fexp, &
                                     varname="FEXP", wformat=wformat)
 
            ! read: Mean topographic index
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, JULES50_struc(n)%jules50%ti_mean, &
                                     varname="TI_MEAN", wformat=wformat)
 
            ! read: Standard deviation of topographic index
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, JULES50_struc(n)%jules50%ti_sig, &
                                     varname="TI_SIG", wformat=wformat)
 
            ! read: Water table depth
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, JULES50_struc(n)%jules50%zw, &
                                     varname="ZW", wformat=wformat)
 
            ! read: soil moist fraction in deep (water table) layer
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, JULES50_struc(n)%jules50%sthzw, &
                                     varname="STHZW", wformat=wformat)
            
            ! read: Wetland fraction
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, JULES50_struc(n)%jules50%fwetl, &
                                     varname="FWETL", wformat=wformat)
            ! read: Surface saturation fraction
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, JULES50_struc(n)%jules50%fsat, &
                                     varname="FSAT", wformat=wformat)
 
            ! read: Unfrozen soil moisture content of the layers as a fraction of saturation
            do l=1, JULES50_struc(n)%sm_levels ! TODO: check loop
                call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, tmptilen, varname="STHU", &
                                         dim=l, vlevels = JULES50_struc(n)%sm_levels, wformat=wformat)
                do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
                    JULES50_struc(n)%jules50(t)%p_s_sthu(l) = tmptilen(t)
                enddo
            enddo
 
            ! read: Frozen soil moisture content of the layers as a fraction of saturation
            do l=1, JULES50_struc(n)%sm_levels ! TODO: check loop
                call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, tmptilen, varname="STHF", &
                                         dim=l, vlevels = JULES50_struc(n)%sm_levels, wformat=wformat)
                do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
                    JULES50_struc(n)%jules50(t)%p_s_sthf(l) = tmptilen(t)
                enddo
            enddo
 
            ! read: Snow layer ice mass on tiles
            do l=1, jules50_struc(n)%temp_tiles*jules50_struc(n)%temp_layers ! TODO: check loop
                call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, tmptilen, varname="SICE", &
                                         dim=l, vlevels = jules50_struc(n)%temp_tiles*jules50_struc(n)%temp_layers, wformat=wformat)
                call l_to_rc(l, jules50_struc(n)%temp_tiles, r_idx, c_idx)
                do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
                    JULES50_struc(n)%jules50(t)%sice(r_idx, c_idx) = tmptilen(t)
                enddo
            enddo
 
            ! read: Snow layer liquid mass on tiles
            do l=1, jules50_struc(n)%temp_tiles*jules50_struc(n)%temp_layers ! TODO: check loop
                call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, tmptilen, varname="SLIQ", &
                                         dim=l, vlevels = jules50_struc(n)%temp_tiles*jules50_struc(n)%temp_layers, wformat=wformat)
                call l_to_rc(l, jules50_struc(n)%temp_tiles, r_idx, c_idx)
                do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
                    JULES50_struc(n)%jules50(t)%sliq(r_idx, c_idx) = tmptilen(t)
                enddo
            enddo
 
            ! read: Snow depth on ground on tiles
            do l=1, JULES50_struc(n)%nsurft ! TODO: check loop
                call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, tmptilen, varname="SNOWDEPTH", &
                                         dim=l, vlevels = JULES50_struc(n)%nsurft, wformat=wformat)
                do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
                    JULES50_struc(n)%jules50(t)%snowdepth(l) = tmptilen(t)
                enddo
            enddo
 
            ! read: Snow layer temperature
            do l=1, jules50_struc(n)%temp_tiles*jules50_struc(n)%temp_layers ! TODO: check loop
                call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, tmptilen, varname="TSNOW", &
                                         dim=l, vlevels = jules50_struc(n)%temp_tiles*jules50_struc(n)%temp_layers, wformat=wformat)
                call l_to_rc(l, jules50_struc(n)%temp_tiles, r_idx, c_idx)
                do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
                    JULES50_struc(n)%jules50(t)%tsnow(r_idx, c_idx) = tmptilen(t)
                enddo
            enddo
 
            ! read: Snow layer grain size on tiles
            do l=1, jules50_struc(n)%temp_tiles*jules50_struc(n)%temp_layers ! TODO: check loop
                call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, tmptilen, varname="RGRAINL", &
                                         dim=l, vlevels = jules50_struc(n)%temp_tiles*jules50_struc(n)%temp_layers, wformat=wformat)
                call l_to_rc(l, jules50_struc(n)%temp_tiles, r_idx, c_idx)
                do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
                    JULES50_struc(n)%jules50(t)%rgrainl(r_idx, c_idx) = tmptilen(t)
                enddo
            enddo
 
            ! read: Snowpack bulk density
            do l=1, JULES50_struc(n)%nsurft ! TODO: check loop
                call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, tmptilen, varname="RHO_SNOW_GRND", &
                                         dim=l, vlevels = JULES50_struc(n)%nsurft, wformat=wformat)
                do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
                    JULES50_struc(n)%jules50(t)%rho_snow_grnd(l) = tmptilen(t)
                enddo
            enddo
 
            ! read: Snow layer densities
            do l=1, JULES50_struc(n)%nsurft*jules50_struc(n)%nsmax ! TODO: check loop
                call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, tmptilen, varname="RHO_SNOW", &
                                         dim=l, vlevels = JULES50_struc(n)%nsurft*JULES50_struc(n)%nsmax, wformat=wformat)
                call l_to_rc(l, jules50_struc(n)%nsurft, r_idx, c_idx)
                do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
                    JULES50_struc(n)%jules50(t)%rho_snow(r_idx, c_idx) = tmptilen(t)
                enddo
            enddo
 
            ! read: Snow layer thickness
            do l=1, jules50_struc(n)%temp_tiles*jules50_struc(n)%temp_layers ! TODO: check loop
                call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, tmptilen, varname="DS", &
                                         dim=l, vlevels = jules50_struc(n)%temp_tiles*jules50_struc(n)%temp_layers, wformat=wformat)
                call l_to_rc(l, jules50_struc(n)%temp_tiles, r_idx, c_idx)
                do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
                    JULES50_struc(n)%jules50(t)%ds(r_idx, c_idx) = tmptilen(t)
                enddo
            enddo

            ! newly added  
           ! Number of snow layers on ground on tiles
            do l=1, JULES50_struc(n)%nsurft ! TODO: check loop
                call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, tmptilen, varname="NSNOW", &
                                         dim=l, vlevels = JULES50_struc(n)%nsurft, wformat=wformat)
                do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
                    JULES50_struc(n)%jules50(t)%nsnow(l) = nint(tmptilen(t))
                enddo
            enddo
            
            ! Tiled land-ice bedrock subsurface temperatures (K)
            do l=1, JULES50_struc(n)%nsurft ! TODO: check loop
                call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, tmptilen, varname="TSURF_ELEV_SURFT", &
                                         dim=l, vlevels = JULES50_struc(n)%nsurft, wformat=wformat)
                do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
                    JULES50_struc(n)%jules50(t)%tsurf_elev_surft(l) = tmptilen(t)
                enddo
            enddo
            
            ! read: Pasture fraction from previous TRIFFID call
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, JULES50_struc(n)%jules50%frac_past_prev, &
                                     varname="FRAC_PAST_PREV", wformat=wformat)
            
            ! read: z0msea
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, JULES50_struc(n)%jules50%z0msea, &
                                     varname="Z0MSEA", wformat=wformat)
            
            ! Gridbox Inorganic N pool on soil levels
            do l=1, JULES50_struc(n)%dim_cslayer ! TODO: check loop
                call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, tmptilen, varname="N_INORG_SOILT_LYRS", &
                                         dim=l, vlevels = JULES50_struc(n)%dim_cslayer, wformat=wformat)
                do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
                    JULES50_struc(n)%jules50(t)%n_inorg_soilt_lyrs(l) = tmptilen(t)
                enddo
            enddo
            
            ! read: Availabile inorganic N for PFTs (depends on roots)  
            do l=1, JULES50_struc(n)%npft*JULES50_struc(n)%dim_cslayer 
                call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, tmptilen, varname="N_INORG_AVAIL_PFT", &
                                         dim=l, vlevels = JULES50_struc(n)%npft*JULES50_struc(n)%dim_cslayer, wformat=wformat)
                call l_to_rc(l, jules50_struc(n)%npft, r_idx, c_idx)
                do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
                    JULES50_struc(n)%jules50(t)%n_inorg_avail_pft(r_idx, c_idx) = tmptilen(t)
                enddo
            enddo
            
            ! read: Atmospheric CO2 fluxes from TRIFFID (kgC/m2/yr)
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, JULES50_struc(n)%jules50%triffid_co2_gb, &
                                     varname="TRIFFID_CO2_GB", wformat=wformat)
            
            ! read: Gridbox canopy water content
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, JULES50_struc(n)%jules50%canopy_gb, &
                                     varname="CANOPY_GB", wformat=wformat)
            
            ! read: Gridbox snowmass (kg/m2)
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, JULES50_struc(n)%jules50%snow_mass_ij, &
                                     varname="SNOW_MASS_IJ", wformat=wformat)
       
            ! close restart file
            if(wformat .eq. "binary") then
                call LIS_releaseUnitNumber(ftn)
            elseif(wformat .eq. "netcdf") then
#if (defined USE_NETCDF3 || defined USE_NETCDF4)
                status = nf90_close(ftn)
                call LIS_verify(status, "Error in nf90_close in jules50_readrst")
#endif
            endif
            deallocate(tmptilen)
        endif    
    enddo
end subroutine jules50_readrst
        
