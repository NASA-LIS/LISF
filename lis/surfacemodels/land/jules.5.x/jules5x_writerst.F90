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
! !ROUTINE: jules5x_writerst
! \label{jules5x_writerst}
! !REVISION HISTORY:
! 16 May 2016; Shugong Wang; initial implementation for JULES 4.3
! 01 Feb 2018; Shugong Wang; updated for JULES 5.0 
! 28 Nov 2018; Shugong Wang; updated for JULES 5.2 
! 12 Dec 2018; Shugong Wang; updated for JULES 5.3
! 08 Jul 2019; Shugong Wang; updated for JULES 5.5
!
! !INTERFACE:
subroutine jules5x_writerst(n)
! !USES:
    use LIS_coreMod, only    : LIS_rc, LIS_masterproc
    use LIS_timeMgrMod, only : LIS_isAlarmRinging
    use LIS_logMod, only     : LIS_logunit, LIS_getNextUnitNumber, &
                               LIS_releaseUnitNumber , LIS_verify
    use LIS_fileIOMod, only  : LIS_create_output_directory, &
                               LIS_create_restart_filename
    use jules5x_lsmMod

#if (defined USE_NETCDF3 || defined USE_NETCDF4)
    use netcdf
#endif

    implicit none
    ! !ARGUMENTS:
    integer, intent(in) :: n
!
! !DESCRIPTION:
!  This program writes restart files for jules5x.
!  This includes all relevant water/energy storage and tile information.
!
!  The routines invoked are:
! \begin{description}
! \item[LIS\_create\_output\_directory](\ref{LIS_create_output_directory})\\
!  creates a timestamped directory for the restart files
! \item[LIS\_create\_restart\_filename](\ref{LIS_create_restart_filename})\\
!  generates a timestamped restart filename
! \item[jules5x\_dump\_restart](\ref{jules5x_dump_restart})\\
!   writes the jules5x variables into the restart file
! \end{description}
!EOP

    character*100 :: filen
    character*20  :: wformat
    logical       :: alarmCheck
    integer       :: ftn, k, l, t, c  
    integer       :: status
    
    ! set restart alarm
    alarmCheck = LIS_isAlarmRinging(LIS_rc, "JULES restart alarm")
    
    ! set restart file format (read from LIS configration file_
    wformat = trim(JULES5x_struc(n)%rformat)
    
    if(alarmCheck .or. (LIS_rc%endtime ==1)) then
        If (LIS_masterproc) Then
            call LIS_create_output_directory("SURFACEMODEL")
            call LIS_create_restart_filename(n, filen, "SURFACEMODEL", "JULES50",&
                                             wformat=wformat)
            if(wformat .eq. "binary") then
                ftn = LIS_getNextUnitNumber()
                open(ftn,file=filen,status="unknown", form="unformatted")
            elseif(wformat .eq. "netcdf") then
#if (defined USE_NETCDF4)
                status = nf90_create(path=filen, cmode=nf90_hdf5, ncid = ftn)
                call LIS_verify(status,"Error in nf90_open in jules5x_writerst")
#endif
#if (defined USE_NETCDF3)
                status = nf90_create(Path = filen, cmode = nf90_clobber, ncid = ftn)
                call LIS_verify(status, "Error in nf90_open in jules5x_writerst")
#endif
             endif
        endif
    
        call jules5x_dump_restart(n, ftn, wformat)
    
        if (LIS_masterproc) then
            if(wformat .eq. "binary") then
                call LIS_releaseUnitNumber(ftn)
            elseif(wformat .eq. "netcdf") then
#if (defined USE_NETCDF3 || defined USE_NETCDF4)
                status = nf90_close(ftn)
                call LIS_verify(status, "Error in nf90_close in jules5x_writerst")
#endif
            endif
            write(LIS_logunit, *) "jules5x archive restart written: ", filen
        endif
    endif
end subroutine jules5x_writerst

!BOP
!
! !ROUTINE: jules5x_dump_restart
! \label{jules5x_dump_restart}
!
! !REVISION HISTORY:
!  This subroutine is generated with the Model Implementation Toolkit developed
!  by Shugong Wang for the NASA Land Information System Version 7. The initial 
!  specification of the subroutine is defined by Sujay Kumar. 
!  5/9/16: Shugong Wang, initial implementation for LIS 7 and jules5x
! !INTERFACE:
subroutine jules5x_dump_restart(n, ftn, wformat)

! !USES:
    use LIS_coreMod, only : LIS_rc, LIS_masterproc
    use LIS_logMod, only  : LIS_logunit
    use LIS_historyMod
    use jules5x_lsmMod

    implicit none

    integer, intent(in) :: ftn
    integer, intent(in) :: n
    character(len=*), intent(in) :: wformat
!
! !DESCRIPTION:
!  This routine gathers the necessary restart variables and performs
!  the actual write statements to create the restart files.
!
!  The arguments are:
!  \begin{description}
!   \item[n]
!    index of the nest
!   \item[ftn]
!    unit number for the restart file
!   \item[wformat]
!    restart file format (binary/netcdf)
!  \end{description}
!
!
!  The following is the list of variables written in the jules5x
!  restart file:
!  \begin{verbatim}
!    nc, nr, nsurfts             - grid and tile space dimensions
!    tsoil_deep                 - jules5x deep soil temperatures [K]
!    wood_prod_fast             - jules5x Fast-turnover wood product C pool [-]
!    wood_prod_med              - jules5x Medium-turnover wood product C pool [-]
!    wood_prod_slow             - jules5x Slow-turnover wood product C pool [-]
!    frac_agr_prev              - jules5x Agricultural fraction from previous TRIFFID call [-]
!    frac_agr                   - jules5x Agricultural fraction [-]
!    n_inorg                    - jules5x Gridbox Inorganic N pool [kg m-2]
!    canht_ft                   - jules5x Canopy height [m]
!    canopy                     - jules5x Surface/canopy water for snow-free land tiles [kg m-2]
!    ns                         - jules5x Soil Organic Nitrogen [kg m-2]
!    cs                         - jules5x Soil carbo [kg m-2]
!    gc                         - jules5x Stomatal conductance to evaporation for land tiles [m s-1]
!    gs                         - jules5x Stomatal conductance to evaporation [m s-1]
!    lai                        - jules5x LAI of plant functional types [-]
!    rgrain                     - jules5x Snow surface grain size on tiles [micron]
!    smc                        - jules5x Soil moisture in a layer at the surface [kg m-2]
!    smcl                       - jules5x Soil moisture content of layers [kg m-2]
!    snow_tile                  - jules5x Lying snow on tiles [kg m-2]
!    snow_grnd                  - jules5x Snow on the ground [kg m-2]
!    soot                       - jules5x Snow soot content (kg kg-1) [kg/kg]
!    t_soil                     - jules5x Sub-surface temperatures [K]
!    tstar_tile                 - jules5x Tile surface temperatures [K]
!    asteps_since_triffid       - jules5x Number of atmospheric timesteps since last call to TRIFFID [-]
!    lai_phen                   - jules5x LAI of PFTs after phenology. Required as separate variable for top-level argument list matching with VEG_IC2A [-]
!    c_veg                      - jules5x Total carbon content of the vegetation [kg m-2]
!    cv                         - jules5x Gridbox mean vegetation carbon [kg m-2]
!    fexp                       - jules5x Decay factor in Sat. Conductivity in water table layer [-]
!    ti_mean                    - jules5x Mean topographic index [-]
!    ti_sig                     - jules5x Standard deviation of topographic index [-]
!    zw                         - jules5x Water table depth [m]
!    sthzw                      - jules5x soil moist fraction in deep (water table) layer []
!    sthu                       - jules5x Unfrozen soil moisture content of the layers as a fraction of saturation [-]
!    sthf                       - jules5x Frozen soil moisture content of the layers as a fraction of saturation [-]
!    sice                       - jules5x Snow layer ice mass on tiles [kg m-2]
!    sliq                       - jules5x Snow layer liquid mass on tiles [kg m-2]
!    snowdepth                  - jules5x Snow depth on ground on tiles [m]
!    tsnow                      - jules5x Snow layer temperature [K]
!    rgrainl                    - jules5x Snow layer grain size on tiles [microns]
!    rho_snow_grnd              - jules5x Snowpack bulk density [kg/m3]
!    rho_snow                   - jules5x Snow layer densities [m]
!    ds                         - jules5x Snow layer thickness [m]
!  \end{verbatim}
!
! The routines invoked are:
! \begin{description}
!   \item[LIS\_writeGlobalHeader\_restart](\ref{LIS_writeGlobalHeader_restart})\\
!      writes the global header information
!   \item[LIS\_writeHeader\_restart](\ref{LIS_writeHeader_restart})\\
!      writes the header information for a variable
!   \item[LIS\_closeHeader\_restart](\ref{LIS_closeHeader_restart})\\
!      close the header
!   \item[LIS\_writevar\_restart](\ref{LIS_writevar_restart})\\
!      writes a variable to the restart file
! \end{description}
! 
!EOP 
               
    integer :: l, t 
    real    :: tmptilen(LIS_rc%npatch(n, LIS_rc%lsm_index))
    integer :: dimID(11)
    integer :: tsoil_deep_ID
    integer :: wood_prod_fast_ID
    integer :: wood_prod_med_ID
    integer :: wood_prod_slow_ID
    integer :: frac_agr_prev_ID
    integer :: frac_past_prev_ID
    integer :: frac_agr_ID
    integer :: n_inorg_ID
    integer :: canht_ft_ID
    integer :: canopy_ID
    integer :: canopy_gb_ID
    integer :: snow_mass_ij_ID
    integer :: ns_ID
    integer :: cs_ID
    integer :: gc_ID
    integer :: gs_ID
    integer :: lai_ID
    integer :: rgrain_ID
    integer :: smc_ID
    integer :: smcl_ID
    integer :: snow_tile_ID
    integer :: snow_grnd_ID
    integer :: soot_ID
    integer :: t_soil_ID
    integer :: tstar_tile_ID
    integer :: asteps_since_triffid_ID
    integer :: lai_phen_ID
    integer :: c_veg_ID
    integer :: cv_ID
    integer :: fexp_ID
    integer :: ti_mean_ID
    integer :: ti_sig_ID
    integer :: zw_ID
    integer :: sthzw_ID
    integer :: fsat_ID 
    integer :: fwetl_ID 
    integer :: sthu_ID
    integer :: sthf_ID
    integer :: sice_ID
    integer :: sliq_ID
    integer :: snowdepth_ID
    integer :: tsnow_ID
    integer :: rgrainl_ID
    integer :: rho_snow_grnd_ID
    integer :: rho_snow_ID
    integer :: ds_ID
    integer :: nsnow_ID 
    integer :: z0msea_ID 
    integer :: tsurf_elev_surft_ID
    integer :: n_inorg_soilt_lyrs_ID 
    integer :: n_inorg_avail_pft_ID 
    integer :: triffid_co2_gb_ID 
    integer :: r_idx, c_idx 
    
    ! write the header of the restart file
    call LIS_writeGlobalHeader_restart(ftn, n, LIS_rc%lsm_index,                       &
                                   "JULES50",                                          & 
                                   dim1 =jules5x_struc(n)%nsurft,                        &   
                                   dim2 =jules5x_struc(n)%npft,                         &
                                   dim3 =jules5x_struc(n)%sm_levels,                        &
                                   dim4 =jules5x_struc(n)%nsmax,                        &
                                   dim5 =jules5x_struc(n)%ns_deep,                      &
                                   dim6 =jules5x_struc(n)%dim_cs1,                      &
                                   dim7 =jules5x_struc(n)%nsurft*jules5x_struc(n)%nsmax, &
                                   dim8 =jules5x_struc(n)%temp_tiles*jules5x_struc(n)%temp_layers, &
                                   dim9 =jules5x_struc(n)%dim_cslayer,                  &
                                   dim10=jules5x_struc(n)%npft*jules5x_struc(n)%dim_cslayer,  &
                                   dimID=dimID,                                        &
                                   output_format = trim(wformat))

    ! write the header for state variable nsnow_surft
    call LIS_writeHeader_restart(ftn, n, dimID, nsnow_ID, "NSNOW",    &
                                 "Number of snow layers on ground on tiles",      &
                                 "-", vlevels=JULES5x_struc(n)%nsurft ,           &
                                 valid_min=-99999.0, valid_max=99999.0,           &
                                 var_flag = "dim1")
    ! Number of snow layers on ground on tiles
    do l=1, JULES5x_struc(n)%nsurft   
        tmptilen = 0
        do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
            tmptilen(t) = 1.0 * JULES5x_struc(n)%jules5x(t)%nsnow(l)
        enddo
        call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, tmptilen, &
                                  varid=nsnow_ID, dim=l, wformat=wformat)
    enddo
    
    ! write the header for state variable tsurf_elev_surft
    call LIS_writeHeader_restart(ftn, n, dimID, tsurf_elev_surft_ID, "TSURF_ELEV_SURFT",    &
                                 "Tiled land-ice bedrock subsurface temperatures",      &
                                 "K", vlevels=JULES5x_struc(n)%nsurft ,           &
                                 valid_min=-99999.0, valid_max=99999.0,           &
                                 var_flag = "dim1")
    ! Number of snow layers on ground on tiles
    do l=1, JULES5x_struc(n)%nsurft   
        tmptilen = 0
        do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
            tmptilen(t) = JULES5x_struc(n)%jules5x(t)%tsurf_elev_surft(l)
        enddo
        call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, tmptilen, &
                                  varid=tsurf_elev_surft_ID, dim=l, wformat=wformat)
    enddo
    ! write the header for state variable z0msea
    call LIS_writeHeader_restart(ftn, n, dimID, z0msea_ID, "Z0MSEA", &
                                 "Sea-surface roughness length for momentum", &
                                 "m", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
    !
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, JULES5x_struc(n)%jules5x%z0msea, &
                              varid=z0msea_ID, dim=1, wformat=wformat)
    
    ! write the header for state variable frac_past_prev
    call LIS_writeHeader_restart(ftn, n, dimID, frac_past_prev_ID, "FRAC_PAST_PREV", &
                                 "Pasture fraction from previous TRIFFID call", &
                                 "-", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
    !Pasture fraction from previous TRIFFID call 
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, JULES5x_struc(n)%jules5x%frac_past_prev, &
                              varid=frac_past_prev_ID, dim=1, wformat=wformat)
    
    ! write the header for state variable n_inorg_soilt_lyrs
    call LIS_writeHeader_restart(ftn, n, dimID, n_inorg_soilt_lyrs_ID, "N_INORG_SOILT_LYRS",    &
                                 "Gridbox Inorganic N pool on soil levels",        &
                                 "kg N/m2", vlevels=JULES5x_struc(n)%dim_cslayer , &
                                 valid_min=-99999.0, valid_max=99999.0,            &
                                 var_flag = "dim9")
    !  Gridbox Inorganic N pool on soil levels (kg N/m2)
    do l=1, JULES5x_struc(n)%dim_cslayer   
        tmptilen = 0
        do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
            tmptilen(t) = 1.0 * JULES5x_struc(n)%jules5x(t)%n_inorg_soilt_lyrs(l)
        enddo
        call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, tmptilen,                 &
                                  varid=n_inorg_soilt_lyrs_ID, dim=l, wformat=wformat)
    enddo
    
    ! write the header for state variable n_inorg_avail_pft
    call LIS_writeHeader_restart(ftn, n, dimID, n_inorg_avail_pft_ID, "N_INORG_AVAIL_PFT",       &
                                 "Availabile inorganic N for PFTs (depends on roots)", "kg N/m2",& 
                                 vlevels=JULES5x_struc(n)%npft*JULES5x_struc(n)%dim_cslayer ,    &
                                 valid_min=-99999.0, valid_max=99999.0,                          &
                                 var_flag = "dim10")
    ! Availabile inorganic N for PFTs (depends on roots) 
    do l=1, JULES5x_struc(n)%npft*JULES5x_struc(n)%dim_cslayer
        call l_to_rc(l, jules5x_struc(n)%npft, r_idx, c_idx)
        tmptilen = 0
        do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
            tmptilen(t) = JULES5x_struc(n)%jules5x(t)%n_inorg_avail_pft(r_idx, c_idx)
        enddo
        call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, tmptilen, &
                                  varid=n_inorg_avail_pft_ID, dim=l, wformat=wformat)
    enddo
   
    ! write the header for state variable triffid_co2_gb
    call LIS_writeHeader_restart(ftn, n, dimID, triffid_co2_gb_ID, "TRIFFID_CO2_GB", &
                                "Atmospheric CO2 fluxes from TRIFFID", &
                                "(kgC/m2/yr)", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
    !Pasture fraction from previous TRIFFID call 
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, JULES5x_struc(n)%jules5x%triffid_co2_gb, &
                             varid=triffid_co2_gb_ID, dim=1, wformat=wformat)
   

    ! write the header for state variable canopy_gb
    call LIS_writeHeader_restart(ftn, n, dimID, canopy_gb_ID, "CANOPY_GB", &
                                "Gridbox canopy water content", &
                                "(kg/m2)", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
    !Pasture fraction from previous TRIFFID call 
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, JULES5x_struc(n)%jules5x%canopy_gb, &
                             varid=canopy_gb_ID, dim=1, wformat=wformat)
    
    ! write the header for state variable snow_mass_ij
    call LIS_writeHeader_restart(ftn, n, dimID, snow_mass_ij_ID, "SNOW_MASS_IJ", &
                                "Gridbox snowmass", &
                                "(kg/m2)", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
    !Pasture fraction from previous TRIFFID call 
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, JULES5x_struc(n)%jules5x%snow_mass_ij, &
                             varid=snow_mass_ij_ID, dim=1, wformat=wformat)
    

    ! write the header for state variable tsoil_deep
    !TODO: check dimension of the state variable following "vlevels="
    !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
    call LIS_writeHeader_restart(ftn, n, dimID, tsoil_deep_ID, "TSOIL_DEEP", &
                                 "deep soil temperatures", &
                                 "K", vlevels=JULES5x_struc(n)%ns_deep , valid_min=-99999.0, valid_max=99999.0, &
                                 var_flag = "dim5") 
 
    ! write the header for state variable wood_prod_fast
    !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
    call LIS_writeHeader_restart(ftn, n, dimID, wood_prod_fast_ID, "WOOD_PROD_FAST", &
                                 "Fast-turnover wood product C pool", &
                                 "-", vlevels=1, valid_min=-99999.0, valid_max=99999.0)

    ! write the header for state variable wood_prod_med
    !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
    call LIS_writeHeader_restart(ftn, n, dimID, wood_prod_med_ID, "WOOD_PROD_MED", &
                                 "Medium-turnover wood product C pool", &
                                 "-", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
    
    ! write the header for state variable wood_prod_slow
    !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
    call LIS_writeHeader_restart(ftn, n, dimID, wood_prod_slow_ID, "WOOD_PROD_SLOW", &
                                 "Slow-turnover wood product C pool", &
                                 "-", vlevels=1, valid_min=-99999.0, valid_max=99999.0)

    ! write the header for state variable frac_agr_prev
    !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
    call LIS_writeHeader_restart(ftn, n, dimID, frac_agr_prev_ID, "FRAC_AGR_PREV", &
                                 "Agricultural fraction from previous TRIFFID call", &
                                 "-", vlevels=1, valid_min=-99999.0, valid_max=99999.0)

    ! write the header for state variable frac_agr
    !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
    call LIS_writeHeader_restart(ftn, n, dimID, frac_agr_ID, "FRAC_AGR", &
                                 "Agricultural fraction", &
                                 "-", vlevels=1, valid_min=-99999.0, valid_max=99999.0)

    ! write the header for state variable n_inorg
    !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
    call LIS_writeHeader_restart(ftn, n, dimID, n_inorg_ID, "N_INORG", &
                                 "Gridbox Inorganic N pool", &
                                 "kg m-2", vlevels=1, valid_min=-99999.0, valid_max=99999.0)

    ! write the header for state variable canht_ft
    !TODO: check dimension of the state variable following "vlevels="
    !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
    call LIS_writeHeader_restart(ftn, n, dimID, canht_ft_ID, "CANHT_FT", &
                                 "Canopy height", &
                                 "m", vlevels=JULES5x_struc(n)%npft , valid_min=-99999.0, valid_max=99999.0, &
                                 var_flag = "dim2")
 
    ! write the header for state variable canopy
    !TODO: check dimension of the state variable following "vlevels="
    !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
    call LIS_writeHeader_restart(ftn, n, dimID, canopy_ID, "CANOPY", &
                                 "Surface/canopy water for snow-free land tiles", &
                                 "kg m-2", vlevels=JULES5x_struc(n)%nsurft , valid_min=-99999.0, valid_max=99999.0, &
                                 var_flag = "dim1")
 
    ! write the header for state variable ns
    !TODO: check dimension of the state variable following "vlevels="
    !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
    call LIS_writeHeader_restart(ftn, n, dimID, ns_ID, "NS", &
                                 "Soil Organic Nitrogen", &
                                 "kg m-2", vlevels=JULES5x_struc(n)%dim_cs1 , valid_min=-99999.0, valid_max=99999.0, &
                                 var_flag = "dim6") 
 
    ! write the header for state variable cs
    !TODO: check dimension of the state variable following "vlevels="
    !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
    call LIS_writeHeader_restart(ftn, n, dimID, cs_ID, "CS", &
                                 "Soil carbo", &
                                 "kg m-2", vlevels=JULES5x_struc(n)%dim_cs1 , valid_min=-99999.0, valid_max=99999.0, &
                                 var_flag = "dim6") 
 
    ! write the header for state variable gc
    !TODO: check dimension of the state variable following "vlevels="
    !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
    call LIS_writeHeader_restart(ftn, n, dimID, gc_ID, "GC", &
                                 "Stomatal conductance to evaporation for land tiles", &
                                 "m s-1", vlevels=JULES5x_struc(n)%nsurft , valid_min=-99999.0, valid_max=99999.0, &
                                 var_flag = "dim1") 
 
    ! write the header for state variable gs
    !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
    call LIS_writeHeader_restart(ftn, n, dimID, gs_ID, "GS", &
                                 "Stomatal conductance to evaporation", &
                                 "m s-1", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
    ! write the header for state variable lai
    !TODO: check dimension of the state variable following "vlevels="
    !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
    call LIS_writeHeader_restart(ftn, n, dimID, lai_ID, "LAI", &
                                 "LAI of plant functional types", &
                                 "-", vlevels=JULES5x_struc(n)%npft , valid_min=-99999.0, valid_max=99999.0, &
                                 var_flag = "dim2") 
 
    ! write the header for state variable rgrain
    !TODO: check dimension of the state variable following "vlevels="
    !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
    call LIS_writeHeader_restart(ftn, n, dimID, rgrain_ID, "RGRAIN", &
                                 "Snow surface grain size on tiles", &
                                 "micron", vlevels=JULES5x_struc(n)%nsurft , valid_min=-99999.0, valid_max=99999.0, &
                                 var_flag = "dim1")
 
    ! write the header for state variable smc
    !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
    call LIS_writeHeader_restart(ftn, n, dimID, smc_ID, "SMC", &
                                 "Soil moisture in a layer at the surface", &
                                 "kg m-2", vlevels=1, valid_min=-99999.0, valid_max=99999.0)

    ! write the header for state variable smcl
    !TODO: check dimension of the state variable following "vlevels="
    !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
    call LIS_writeHeader_restart(ftn, n, dimID, smcl_ID, "SMCL", &
                                 "Soil moisture content of layers", &
                                 "kg m-2", vlevels=JULES5x_struc(n)%sm_levels , valid_min=-99999.0, valid_max=99999.0, &
                                 var_flag = "dim3") 
 
    ! write the header for state variable snow_tile
    !TODO: check dimension of the state variable following "vlevels="
    !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
    call LIS_writeHeader_restart(ftn, n, dimID, snow_tile_ID, "SNOW_TILE", &
                                 "Lying snow on tiles", &
                                 "kg m-2", vlevels=JULES5x_struc(n)%nsurft , valid_min=-99999.0, valid_max=99999.0, &
                                 var_flag = "dim1") 
 
    ! write the header for state variable snow_grnd
    !TODO: check dimension of the state variable following "vlevels="
    !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
    call LIS_writeHeader_restart(ftn, n, dimID, snow_grnd_ID, "SNOW_GRND", &
                                 "Snow on the ground", &
                                 "kg m-2", vlevels=JULES5x_struc(n)%nsurft , valid_min=-99999.0, valid_max=99999.0, &
                                 var_flag = "dim1") 
 
    ! write the header for state variable soot
    !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
    call LIS_writeHeader_restart(ftn, n, dimID, soot_ID, "SOOT", &
                                 "Snow soot content (kg kg-1)", &
                                 "kg kg-1", vlevels=1, valid_min=-99999.0, valid_max=99999.0)

    ! write the header for state variable t_soil
    !TODO: check dimension of the state variable following "vlevels="
    !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
    call LIS_writeHeader_restart(ftn, n, dimID, t_soil_ID, "T_SOIL", &
                                 "Sub-surface temperatures", &
                                 "K", vlevels=JULES5x_struc(n)%sm_levels , valid_min=-99999.0, valid_max=99999.0, &
                                 var_flag = "dim3")
 
    ! write the header for state variable tstar_tile
    !TODO: check dimension of the state variable following "vlevels="
    !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
    call LIS_writeHeader_restart(ftn, n, dimID, tstar_tile_ID, "TSTAR_TILE", &
                                 "Tile surface temperatures", &
                                 "K", vlevels=JULES5x_struc(n)%nsurft , valid_min=-99999.0, valid_max=99999.0, &
                                 var_flag = "dim1")
 
    ! write the header for state variable asteps_since_triffid
    !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
    call LIS_writeHeader_restart(ftn, n, dimID, asteps_since_triffid_ID, "ASTEPS_SINCE_TRIFFID", &
                                 "Number of atmospheric timesteps since last call to TRIFFID", &
                                 "-", vlevels=1, valid_min=-99999.0, valid_max=99999.0)

    ! write the header for state variable lai_phen
    !TODO: check dimension of the state variable following "vlevels="
    !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
    call LIS_writeHeader_restart(ftn, n, dimID, lai_phen_ID, "LAI_PHEN", &
                                 "LAI of PFTs after phenology. Required as separate variable for top-level argument list matching with VEG_IC2A", &
                                 "-", vlevels=JULES5x_struc(n)%npft , valid_min=-99999.0, valid_max=99999.0, &
                                 var_flag = "dim2") 
 
    ! write the header for state variable c_veg
    !TODO: check dimension of the state variable following "vlevels="
    !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
    call LIS_writeHeader_restart(ftn, n, dimID, c_veg_ID, "C_VEG", &
                                 "Total carbon content of the vegetation", &
                                 "kg m-2", vlevels=JULES5x_struc(n)%npft , valid_min=-99999.0, valid_max=99999.0, &
                                 var_flag = "dim2") 
 
    ! write the header for state variable cv
    !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
    call LIS_writeHeader_restart(ftn, n, dimID, cv_ID, "CV", &
                                 "Gridbox mean vegetation carbon", &
                                 "kg m-2", vlevels=1, valid_min=-99999.0, valid_max=99999.0)

    ! write the header for state variable fexp
    !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
    call LIS_writeHeader_restart(ftn, n, dimID, fexp_ID, "FEXP", &
                                 "Decay factor in Sat. Conductivity in water table layer", &
                                 "-", vlevels=1, valid_min=-99999.0, valid_max=99999.0)

    ! write the header for state variable ti_mean
    !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
    call LIS_writeHeader_restart(ftn, n, dimID, ti_mean_ID, "TI_MEAN", &
                                 "Mean topographic index", &
                                 "-", vlevels=1, valid_min=-99999.0, valid_max=99999.0)

    ! write the header for state variable ti_sig
    !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
    call LIS_writeHeader_restart(ftn, n, dimID, ti_sig_ID, "TI_SIG", &
                                 "Standard deviation of topographic index", &
                                 "-", vlevels=1, valid_min=-99999.0, valid_max=99999.0)

    ! write the header for state variable zw
    !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
    call LIS_writeHeader_restart(ftn, n, dimID, zw_ID, "ZW", &
                                 "Water table depth", &
                                 "m", vlevels=1, valid_min=-99999.0, valid_max=99999.0)

    ! write the header for state variable sthzw
    !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
    call LIS_writeHeader_restart(ftn, n, dimID, sthzw_ID, "STHZW", &
                                 "soil moist fraction in deep (water table) layer", &
                                 "-", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
    
    ! write the header for state variable fwetl 
    !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
    call LIS_writeHeader_restart(ftn, n, dimID, fwetl_ID, "FWETL", &
                                 "Wetland fraction", &
                                 "-", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
    
    ! write the header for state variable fsat
    !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
    call LIS_writeHeader_restart(ftn, n, dimID, fsat_ID, "FSAT", &
                                 "Surface saturation fraction", &
                                 "-", vlevels=1, valid_min=-99999.0, valid_max=99999.0)

    ! write the header for state variable sthu
    !TODO: check dimension of the state variable following "vlevels="
    !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
    call LIS_writeHeader_restart(ftn, n, dimID, sthu_ID, "STHU", &
                                 "Unfrozen soil moisture content of the layers as a fraction of saturation", &
                                 "-", vlevels=JULES5x_struc(n)%sm_levels , valid_min=-99999.0, valid_max=99999.0, &
                                 var_flag = "dim3")
 
    ! write the header for state variable sthf
    !TODO: check dimension of the state variable following "vlevels="
    !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
    call LIS_writeHeader_restart(ftn, n, dimID, sthf_ID, "STHF", &
                                 "Frozen soil moisture content of the layers as a fraction of saturation", &
                                 "-", vlevels=JULES5x_struc(n)%sm_levels , valid_min=-99999.0, valid_max=99999.0, &
                                 var_flag = "dim3") 
 
    ! write the header for state variable sice
    !TODO: check dimension of the state variable following "vlevels="
    !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
    call LIS_writeHeader_restart(ftn, n, dimID, sice_ID, "SICE", &
                                 "Snow layer ice mass on tiles", &
                                 !"kg m-2", vlevels=JULES5x_struc(n)%nsurft*JULES5x_struc(n)%nsmax , valid_min=-99999.0, valid_max=99999.0, &
                                 "kg m-2", vlevels=jules5x_struc(n)%temp_tiles*jules5x_struc(n)%temp_layers , valid_min=-99999.0, valid_max=99999.0, &
                                 var_flag = "dim8")
 
    ! write the header for state variable sliq
    !TODO: check dimension of the state variable following "vlevels="
    !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
    call LIS_writeHeader_restart(ftn, n, dimID, sliq_ID, "SLIQ", &
                                 "Snow layer liquid mass on tiles", &
                                 !"kg m-2", vlevels=JULES5x_struc(n)%nsurft*JULES5x_struc(n)%nsmax , valid_min=-99999.0, valid_max=99999.0, &
                                 "kg m-2", vlevels=jules5x_struc(n)%temp_tiles*jules5x_struc(n)%temp_layers , valid_min=-99999.0, valid_max=99999.0, &
                                 var_flag = "dim8")
 
    ! write the header for state variable snowdepth
    !TODO: check dimension of the state variable following "vlevels="
    !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
    call LIS_writeHeader_restart(ftn, n, dimID, snowdepth_ID, "SNOWDEPTH", &
                                 "Snow depth on ground on tiles", &
                                 "m", vlevels=JULES5x_struc(n)%nsurft, valid_min=-99999.0, valid_max=99999.0, &
                                 var_flag = "dim1")
 
    ! write the header for state variable tsnow
    !TODO: check dimension of the state variable following "vlevels="
    !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
    call LIS_writeHeader_restart(ftn, n, dimID, tsnow_ID, "TSNOW", &
                                 "Snow layer temperature", &
                                 !"K", vlevels=JULES5x_struc(n)%nsurft*JULES5x_struc(n)%nsmax , valid_min=-99999.0, valid_max=99999.0, &
                                 "K", vlevels=jules5x_struc(n)%temp_tiles*jules5x_struc(n)%temp_layers , valid_min=-99999.0, valid_max=99999.0, &
                                 var_flag = "dim8") 
 
    ! write the header for state variable rgrainl
    !TODO: check dimension of the state variable following "vlevels="
    !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
    call LIS_writeHeader_restart(ftn, n, dimID, rgrainl_ID, "RGRAINL", &
                                 "Snow layer grain size on tiles", &
                                 !"microns", vlevels=JULES5x_struc(n)%nsurft*JULES5x_struc(n)%nsmax , valid_min=-99999.0, valid_max=99999.0, &
                                 "microns", vlevels=jules5x_struc(n)%temp_tiles*jules5x_struc(n)%temp_layers , valid_min=-99999.0, valid_max=99999.0, &
                                 var_flag = "dim8")
 
    ! write the header for state variable rho_snow_grnd
    !TODO: check dimension of the state variable following "vlevels="
    !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
    call LIS_writeHeader_restart(ftn, n, dimID, rho_snow_grnd_ID, "RHO_SNOW_GRND", &
                                 "Snowpack bulk density", &
                                 "kg/m3", vlevels=JULES5x_struc(n)%nsurft , valid_min=-99999.0, valid_max=99999.0, &
                                 var_flag = "dim7")
 
    ! write the header for state variable rho_snow
    !TODO: check dimension of the state variable following "vlevels="
    !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
    call LIS_writeHeader_restart(ftn, n, dimID, rho_snow_ID, "RHO_SNOW", &
                                 "Snow layer densities", &
                                 "m", vlevels=JULES5x_struc(n)%nsurft*JULES5x_struc(n)%nsmax , valid_min=-99999.0, valid_max=99999.0, &
                                 var_flag = "dim7")
 
    ! write the header for state variable ds
    !TODO: check dimension of the state variable following "vlevels="
    !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
    call LIS_writeHeader_restart(ftn, n, dimID, ds_ID, "DS",                                            &
                                 "Snow layer thickness", "m",                                           &
                                 vlevels=jules5x_struc(n)%temp_tiles*jules5x_struc(n)%temp_layers,      & 
                                 valid_min=-99999.0, valid_max=99999.0,                                 &
                                 var_flag = "dim8")
 
    ! close header of restart file
    call LIS_closeHeader_restart(ftn, n, LIS_rc%lsm_index, dimID, JULES5x_struc(n)%rstInterval)

    ! write state variables into restart file
    ! deep soil temperatures
    do l=1, JULES5x_struc(n)%ns_deep   ! TODO: check loop
        tmptilen = 0
        do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
            tmptilen(t) = JULES5x_struc(n)%jules5x(t)%tsoil_deep(l)
        enddo
        call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, tmptilen, &
                                  varid=tsoil_deep_ID, dim=l, wformat=wformat)
    enddo
    ! Fast-turnover wood product C pool
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, JULES5x_struc(n)%jules5x%wood_prod_fast, &
                              varid=wood_prod_fast_ID, dim=1, wformat=wformat)

    ! Medium-turnover wood product C pool
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, JULES5x_struc(n)%jules5x%wood_prod_med, &
                              varid=wood_prod_med_ID, dim=1, wformat=wformat)

    ! Slow-turnover wood product C pool
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, JULES5x_struc(n)%jules5x%wood_prod_slow, &
                              varid=wood_prod_slow_ID, dim=1, wformat=wformat)

    ! Agricultural fraction from previous TRIFFID call
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, JULES5x_struc(n)%jules5x%frac_agr_prev, &
                              varid=frac_agr_prev_ID, dim=1, wformat=wformat)

    ! Agricultural fraction
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, JULES5x_struc(n)%jules5x%frac_agr, &
                              varid=frac_agr_ID, dim=1, wformat=wformat)

    ! Gridbox Inorganic N pool
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, JULES5x_struc(n)%jules5x%n_inorg, &
                              varid=n_inorg_ID, dim=1, wformat=wformat)

    ! Canopy height
    do l=1, JULES5x_struc(n)%npft   ! TODO: check loop
        tmptilen = 0
        do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
            tmptilen(t) = JULES5x_struc(n)%jules5x(t)%canht_ft(l)
        enddo
        call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, tmptilen, &
                                  varid=canht_ft_ID, dim=l, wformat=wformat)
    enddo

    ! Surface/canopy water for snow-free land tiles
    do l=1, JULES5x_struc(n)%nsurft   ! TODO: check loop
        tmptilen = 0
        do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
            tmptilen(t) = JULES5x_struc(n)%jules5x(t)%canopy(l)
        enddo
        call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, tmptilen, &
                                  varid=canopy_ID, dim=l, wformat=wformat)
    enddo
    
    ! Soil Organic Nitrogen
    do l=1, JULES5x_struc(n)%dim_cs1*JULES5x_struc(n)%dim_cslayer   ! TODO: check loop
        call l_to_rc(l, jules5x_struc(n)%dim_cs1, r_idx, c_idx)
        tmptilen = 0
        do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
            tmptilen(t) = JULES5x_struc(n)%jules5x(t)%ns(r_idx,c_idx)
        enddo
        call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, tmptilen, &
                                  varid=ns_ID, dim=l, wformat=wformat)
    enddo
    
    ! Soil carbo
    do l=1, JULES5x_struc(n)%dim_cs1*JULES5x_struc(n)%dim_cslayer   ! TODO: check loop
        call l_to_rc(l, jules5x_struc(n)%dim_cs1, r_idx, c_idx)
        tmptilen = 0
        do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
            tmptilen(t) = JULES5x_struc(n)%jules5x(t)%cs_pool_soilt(r_idx, c_idx)
        enddo
        call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, tmptilen, &
                                  varid=cs_ID, dim=l, wformat=wformat)
    enddo
   
    ! Stomatal conductance to evaporation for land tiles
    do l=1, JULES5x_struc(n)%nsurft   ! TODO: check loop
        tmptilen = 0
        do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
            tmptilen(t) = JULES5x_struc(n)%jules5x(t)%gc_surft(l)
        enddo
        call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, tmptilen, &
                                  varid=gc_ID, dim=l, wformat=wformat)
    enddo
   
    ! Stomatal conductance to evaporation
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, JULES5x_struc(n)%jules5x%gs_gb, &
                              varid=gs_ID, dim=1, wformat=wformat)

   
    ! LAI of plant functional types
    do l=1, JULES5x_struc(n)%npft   ! TODO: check loop
        tmptilen = 0
        do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
            tmptilen(t) = JULES5x_struc(n)%jules5x(t)%lai(l)
        enddo
        call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, tmptilen, &
                                  varid=lai_ID, dim=l, wformat=wformat)
    enddo

    ! Snow surface grain size on tiles
    do l=1, JULES5x_struc(n)%nsurft   ! TODO: check loop
        tmptilen = 0
        do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
            tmptilen(t) = JULES5x_struc(n)%jules5x(t)%rgrain(l)
        enddo
        call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, tmptilen, &
                                  varid=rgrain_ID, dim=l, wformat=wformat)
    enddo

    ! Soil moisture in a layer at the surface
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, JULES5x_struc(n)%jules5x%smc_soilt, &
                              varid=smc_ID, dim=1, wformat=wformat)

    ! Soil moisture content of layers
    do l=1, JULES5x_struc(n)%sm_levels   ! TODO: check loop
        tmptilen = 0
        do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
            tmptilen(t) = JULES5x_struc(n)%jules5x(t)%smcl_soilt(l)
        enddo
        call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, tmptilen, &
                                  varid=smcl_ID, dim=l, wformat=wformat)
    enddo

    ! Lying snow on tiles
    do l=1, JULES5x_struc(n)%nsurft   ! TODO: check loop
        tmptilen = 0
        do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
            tmptilen(t) = JULES5x_struc(n)%jules5x(t)%snow_tile(l)
        enddo
        call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, tmptilen, &
                                  varid=snow_tile_ID, dim=l, wformat=wformat)
    enddo

    ! Snow on the ground
    do l=1, JULES5x_struc(n)%nsurft   ! TODO: check loop
        tmptilen = 0
        do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
            tmptilen(t) = JULES5x_struc(n)%jules5x(t)%snow_grnd(l)
        enddo
        call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, tmptilen, &
                                  varid=snow_grnd_ID, dim=l, wformat=wformat)
    enddo

    ! Snow soot content (kg kg-1)
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, JULES5x_struc(n)%jules5x%soot_ij, &
                              varid=soot_ID, dim=1, wformat=wformat)

    ! Sub-surface temperatures
    do l=1, JULES5x_struc(n)%sm_levels   ! TODO: check loop
        tmptilen = 0
        do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
            tmptilen(t) = JULES5x_struc(n)%jules5x(t)%t_soil(l)
        enddo
        call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, tmptilen, &
                                  varid=t_soil_ID, dim=l, wformat=wformat)
    enddo

    ! Tile surface temperatures
    do l=1, JULES5x_struc(n)%nsurft   ! TODO: check loop
        tmptilen = 0
        do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
            tmptilen(t) = JULES5x_struc(n)%jules5x(t)%tstar_tile(l)
        enddo
        call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, tmptilen, &
                                  varid=tstar_tile_ID, dim=l, wformat=wformat)
    enddo

    ! Number of atmospheric timesteps since last call to TRIFFID
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, JULES5x_struc(n)%jules5x%asteps_since_triffid, &
                              varid=asteps_since_triffid_ID, dim=1, wformat=wformat)

    ! LAI of PFTs after phenology. Required as separate variable for top-level argument list matching with VEG_IC2A
    do l=1, JULES5x_struc(n)%npft   ! TODO: check loop
        tmptilen = 0
        do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
            tmptilen(t) = JULES5x_struc(n)%jules5x(t)%lai_phen(l)
        enddo
        call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, tmptilen, &
                                  varid=lai_phen_ID, dim=l, wformat=wformat)
    enddo

    ! Total carbon content of the vegetation
    do l=1, JULES5x_struc(n)%npft   ! TODO: check loop
        tmptilen = 0
        do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
            tmptilen(t) = JULES5x_struc(n)%jules5x(t)%c_veg(l)
        enddo
        call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, tmptilen, &
                                  varid=c_veg_ID, dim=l, wformat=wformat)
    enddo

    ! Gridbox mean vegetation carbon
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, JULES5x_struc(n)%jules5x%cv, &
                              varid=cv_ID, dim=1, wformat=wformat)

    ! Decay factor in Sat. Conductivity in water table layer
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, JULES5x_struc(n)%jules5x%fexp, &
                              varid=fexp_ID, dim=1, wformat=wformat)

    ! Mean topographic index
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, JULES5x_struc(n)%jules5x%ti_mean, &
                              varid=ti_mean_ID, dim=1, wformat=wformat)

    ! Standard deviation of topographic index
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, JULES5x_struc(n)%jules5x%ti_sig, &
                              varid=ti_sig_ID, dim=1, wformat=wformat)

    ! Water table depth
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, JULES5x_struc(n)%jules5x%zw, &
                              varid=zw_ID, dim=1, wformat=wformat)
    
    ! Surface saturation fraction
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, JULES5x_struc(n)%jules5x%fsat, &
                              varid=fsat_ID, dim=1, wformat=wformat)

    ! Wetland fraction 
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, JULES5x_struc(n)%jules5x%fwetl, &
                              varid=fwetl_ID, dim=1, wformat=wformat)

    ! soil moist fraction in deep (water table) layer
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, JULES5x_struc(n)%jules5x%sthzw, &
                              varid=sthzw_ID, dim=1, wformat=wformat)

    ! Unfrozen soil moisture content of the layers as a fraction of saturation
    do l=1, JULES5x_struc(n)%sm_levels   ! TODO: check loop
        tmptilen = 0
        do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
            tmptilen(t) = JULES5x_struc(n)%jules5x(t)%p_s_sthu(l)
        enddo
        call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, tmptilen, &
                                  varid=sthu_ID, dim=l, wformat=wformat)
    enddo

    ! Frozen soil moisture content of the layers as a fraction of saturation
    do l=1, JULES5x_struc(n)%sm_levels   ! TODO: check loop
        tmptilen = 0
        do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
            tmptilen(t) = JULES5x_struc(n)%jules5x(t)%p_s_sthf(l)
        enddo
        call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, tmptilen, &
                                  varid=sthf_ID, dim=l, wformat=wformat)
    enddo

    ! Snow layer ice mass on tiles
    do l=1, jules5x_struc(n)%temp_tiles*jules5x_struc(n)%temp_layers   ! TODO: check loop
        tmptilen = 0
        call l_to_rc(l, jules5x_struc(n)%temp_tiles, r_idx, c_idx)
        do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
            tmptilen(t) = JULES5x_struc(n)%jules5x(t)%sice(r_idx, c_idx)
        enddo
        call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, tmptilen, &
                                  varid=sice_ID, dim=l, wformat=wformat)
    enddo

    ! Snow layer liquid mass on tiles
    do l=1, jules5x_struc(n)%temp_tiles*jules5x_struc(n)%temp_layers   ! TODO: check loop
        tmptilen = 0
        call l_to_rc(l, jules5x_struc(n)%temp_tiles, r_idx, c_idx)
        do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
            tmptilen(t) = JULES5x_struc(n)%jules5x(t)%sliq(r_idx, c_idx)
        enddo
        call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, tmptilen, &
                                  varid=sliq_ID, dim=l, wformat=wformat)
    enddo

    ! Snow depth on ground on tiles
    do l=1, JULES5x_struc(n)%nsurft   ! TODO: check loop
        tmptilen = 0
        do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
            tmptilen(t) = JULES5x_struc(n)%jules5x(t)%snowdepth(l)
        enddo
        call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, tmptilen, &
                                  varid=snowdepth_ID, dim=l, wformat=wformat)
    enddo

    ! Snow layer temperature
    do l=1, jules5x_struc(n)%temp_tiles*jules5x_struc(n)%temp_layers   ! TODO: check loop
        tmptilen = 0
        call l_to_rc(l, jules5x_struc(n)%temp_tiles, r_idx, c_idx)
        do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
            tmptilen(t) = JULES5x_struc(n)%jules5x(t)%tsnow(r_idx, c_idx)
        enddo
        call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, tmptilen, &
                                  varid=tsnow_ID, dim=l, wformat=wformat)
    enddo

    ! Snow layer grain size on tiles
    do l=1, jules5x_struc(n)%temp_tiles*jules5x_struc(n)%temp_layers   ! TODO: check loop
        tmptilen = 0
        call l_to_rc(l, jules5x_struc(n)%temp_tiles, r_idx, c_idx)
        do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
            tmptilen(t) = JULES5x_struc(n)%jules5x(t)%rgrainl(r_idx, c_idx)
        enddo
        call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, tmptilen, &
                                  varid=rgrainl_ID, dim=l, wformat=wformat)
    enddo

    ! Snowpack bulk density
    do l=1, JULES5x_struc(n)%nsurft   ! TODO: check loop
        tmptilen = 0
        do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
            tmptilen(t) = JULES5x_struc(n)%jules5x(t)%rho_snow_grnd(l)
        enddo
        call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, tmptilen, &
                                  varid=rho_snow_grnd_ID, dim=l, wformat=wformat)
    enddo

    ! Snow layer densities
    do l=1, JULES5x_struc(n)%nsurft*jules5x_struc(n)%nsmax   ! TODO: check loop
        tmptilen = 0
        call l_to_rc(l, jules5x_struc(n)%nsurft, r_idx, c_idx)
        do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
            tmptilen(t) = JULES5x_struc(n)%jules5x(t)%rho_snow(r_idx, c_idx)
        enddo
        call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, tmptilen, &
                                  varid=rho_snow_ID, dim=l, wformat=wformat)
    enddo
   
    ! Snow layer thickness
    do l=1, jules5x_struc(n)%temp_tiles*jules5x_struc(n)%temp_layers   ! TODO: check loop
        call l_to_rc(l, jules5x_struc(n)%temp_tiles, r_idx, c_idx)
        tmptilen = 0
        do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
            tmptilen(t) = JULES5x_struc(n)%jules5x(t)%ds(r_idx, c_idx)
        enddo
        call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, tmptilen, &
                                  varid=ds_ID, dim=l, wformat=wformat)
    enddo
end subroutine jules5x_dump_restart

