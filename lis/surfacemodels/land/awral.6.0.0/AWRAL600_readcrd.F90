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
! !ROUTINE: AWRAL600_readcrd
! \label{AWRAL600\_readcrd}
!
! !REVISION HISTORY:
!  This subroutine is generated with the Model Implementation Toolkit developed
!  by Shugong Wang for the NASA Land Information System Version 7. The initial 
!  specification of the subroutine is defined by Sujay Kumar. 
!   12/18/18 : Wendy Sharples, Shugong Wang, initial implementation for LIS 7 and AWRAL600
!
! !INTERFACE:
subroutine AWRAL600_readcrd()
! !USES:
    use ESMF
    use LIS_coreMod, only    : LIS_rc , LIS_config
    use LIS_timeMgrMod, only : LIS_parseTimeString
    use LIS_logMod, only     : LIS_logunit , LIS_verify
    use AWRAL600_lsmMod, only       : AWRAL600_struc

!
! !DESCRIPTION:
!
!  This routine reads the options specific to AWRAL600 model from
!  the LIS configuration file.
!
!EOP
    implicit none

    integer      :: rc 
    integer      :: n, i
    character*10 :: time 
    character*6  :: str_i
 
    write(LIS_logunit, *) "[INFO] Start reading LIS configuration file for AWRAL600 model"
    
    call ESMF_ConfigFindLabel(LIS_config, "AWRAL600 model timestep:", rc = rc)
    do n=1,LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, Time, rc = rc)
        call LIS_verify(rc, "AWRAL600 model timestep: not defined")
        call LIS_parseTimeString(time, AWRAL600_struc(n)%ts)
    enddo
    
    call ESMF_ConfigFindLabel(LIS_config, "AWRAL600 restart output interval:", rc = rc)
    do n=1,LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, Time, rc = rc)
        call LIS_verify(rc,"AWRAL600 restart output interval: not defined")
        call LIS_parseTimeString(time, AWRAL600_struc(n)%rstInterval)
    enddo
    
    !---------------------------!
    ! Constant Parameters       !
    !---------------------------!
    ! number of hydrologic response units
    call ESMF_ConfigFindLabel(LIS_config, "AWRAL600 nhru:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, AWRAL600_struc(n)%nhru, rc=rc)
        call LIS_verify(rc, "AWRAL600 nhru: not defined")
    enddo

    ! number of hypsometric curve percentile bins
    call ESMF_ConfigFindLabel(LIS_config, "AWRAL600 nhypsbins:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, AWRAL600_struc(n)%nhypsbins, rc=rc)
        call LIS_verify(rc, "AWRAL600 nhypsbins: not defined")
    enddo

    ! scaling factor for slope
    call ESMF_ConfigFindLabel(LIS_config, "AWRAL600 slope_coeff:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, AWRAL600_struc(n)%slope_coeff, rc=rc)
        call LIS_verify(rc, "AWRAL600 slope_coeff: not defined")
    enddo
 
    ! air pressure
    call ESMF_ConfigFindLabel(LIS_config, "AWRAL600 pair:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, AWRAL600_struc(n)%pair, rc=rc)
        call LIS_verify(rc, "AWRAL600 pair: not defined")
    enddo
 
    ! scaling factor for ratio of saturated hydraulic conductivity
    call ESMF_ConfigFindLabel(LIS_config, "AWRAL600 kr_coeff:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, AWRAL600_struc(n)%kr_coeff, rc=rc)
        call LIS_verify(rc, "AWRAL600 kr_coeff: not defined")
    enddo

    ! hypsometric percentile bins
    call ESMF_ConfigFindLabel(LIS_config, "AWRAL600 hypsperc:", rc = rc)
    do n=1, LIS_rc%nnest
        allocate(AWRAL600_struc(n)%hypsperc(AWRAL600_struc(n)%nhypsbins))
        do i = 1, AWRAL600_struc(n)%nhypsbins
            call ESMF_ConfigGetAttribute(LIS_config, AWRAL600_struc(n)%hypsperc(i), rc=rc)
            call LIS_verify(rc, "AWRAL600 hypsperc: not defined")
        enddo
    enddo

    ! dry soil albedo for each hru
    call ESMF_ConfigFindLabel(LIS_config, "AWRAL600 alb_dry:", rc = rc)
    do n=1, LIS_rc%nnest
        allocate(AWRAL600_struc(n)%alb_dry(AWRAL600_struc(n)%nhru))
        do i = 1, AWRAL600_struc(n)%nhru
            call ESMF_ConfigGetAttribute(LIS_config, AWRAL600_struc(n)%alb_dry(i), rc=rc)
            call LIS_verify(rc, "AWRAL600 alb_dry: not defined")
        enddo
    enddo
 
    ! wet soil albedo for each hru
    call ESMF_ConfigFindLabel(LIS_config, "AWRAL600 alb_wet:", rc = rc)
    do n=1, LIS_rc%nnest
        allocate(AWRAL600_struc(n)%alb_wet(AWRAL600_struc(n)%nhru))    
        do i = 1, AWRAL600_struc(n)%nhru
            call ESMF_ConfigGetAttribute(LIS_config, AWRAL600_struc(n)%alb_wet(i), rc=rc)
            call LIS_verify(rc, "AWRAL600 alb_wet: not defined")
        enddo
    enddo
 
    ! coefficient relating vegetation photosynthetic capacity to maximum stomatal conductance for each hru
    call ESMF_ConfigFindLabel(LIS_config, "AWRAL600 cgsmax:", rc = rc)
    do n=1, LIS_rc%nnest
        allocate(AWRAL600_struc(n)%cgsmax(AWRAL600_struc(n)%nhru))
        do i = 1, AWRAL600_struc(n)%nhru
            call ESMF_ConfigGetAttribute(LIS_config, AWRAL600_struc(n)%cgsmax(i), rc=rc)
            call LIS_verify(rc, "AWRAL600 cgsmax: not defined")
        enddo
    enddo
 
    ! specific ratio of the mean evaporation rate and the mean rainfall intensity during storms for each hru
    call ESMF_ConfigFindLabel(LIS_config, "AWRAL600 er_frac_ref:", rc = rc)
    do n=1, LIS_rc%nnest
        allocate(AWRAL600_struc(n)%er_frac_ref(AWRAL600_struc(n)%nhru))
        do i = 1, AWRAL600_struc(n)%nhru
            call ESMF_ConfigGetAttribute(LIS_config, AWRAL600_struc(n)%er_frac_ref(i), rc=rc)
            call LIS_verify(rc, "AWRAL600 er_frac_ref: not defined")
        enddo
    enddo
 
    ! soil evaporation scaling factor corresponding to unlimited soil water supply for each hru
    call ESMF_ConfigFindLabel(LIS_config, "AWRAL600 fsoilemax:", rc = rc)
    do n=1, LIS_rc%nnest
        allocate(AWRAL600_struc(n)%fsoilemax(AWRAL600_struc(n)%nhru))
        do i = 1, AWRAL600_struc(n)%nhru
            call ESMF_ConfigGetAttribute(LIS_config, AWRAL600_struc(n)%fsoilemax(i), rc=rc)
            call LIS_verify(rc, "AWRAL600 fsoilemax: not defined")
        enddo
    enddo
 
    ! reference leaf area index (at which fv = 0.63) for each hru
    call ESMF_ConfigFindLabel(LIS_config, "AWRAL600 lairef:", rc = rc)
    do n=1, LIS_rc%nnest
        allocate(AWRAL600_struc(n)%lairef(AWRAL600_struc(n)%nhru))
        do i = 1, AWRAL600_struc(n)%nhru
            call ESMF_ConfigGetAttribute(LIS_config, AWRAL600_struc(n)%lairef(i), rc=rc)
            call LIS_verify(rc, "AWRAL600 lairef: not defined")
        enddo
    enddo
 
    ! rooting depth for each hru
    call ESMF_ConfigFindLabel(LIS_config, "AWRAL600 rd:", rc = rc)
    do n=1, LIS_rc%nnest
        allocate(AWRAL600_struc(n)%rd(AWRAL600_struc(n)%nhru))
        do i = 1, AWRAL600_struc(n)%nhru
            call ESMF_ConfigGetAttribute(LIS_config, AWRAL600_struc(n)%rd(i), rc=rc)
            call LIS_verify(rc, "AWRAL600 rd: not defined")
        enddo
    enddo
 
    ! specific canopy rainfall storage per unit leaf area for each hru
    call ESMF_ConfigFindLabel(LIS_config, "AWRAL600 s_sls:", rc = rc)
    do n=1, LIS_rc%nnest
        allocate(AWRAL600_struc(n)%s_sls(AWRAL600_struc(n)%nhru))
        do i = 1, AWRAL600_struc(n)%nhru
            call ESMF_ConfigGetAttribute(LIS_config, AWRAL600_struc(n)%s_sls(i), rc=rc)
            call LIS_verify(rc, "AWRAL600 s_sls: not defined")
        enddo
    enddo
 
    ! specific leaf area for each hru
    call ESMF_ConfigFindLabel(LIS_config, "AWRAL600 sla:", rc = rc)
    do n=1, LIS_rc%nnest
        allocate(AWRAL600_struc(n)%sla(AWRAL600_struc(n)%nhru))
        do i = 1, AWRAL600_struc(n)%nhru
            call ESMF_ConfigGetAttribute(LIS_config, AWRAL600_struc(n)%sla(i), rc=rc)
            call LIS_verify(rc, "AWRAL600 sla: not defined")
        enddo
    enddo
 
    ! characteristic time scale for vegetation growth towards equilibrium for each hru
    call ESMF_ConfigFindLabel(LIS_config, "AWRAL600 tgrow:", rc = rc)
    do n=1, LIS_rc%nnest
        allocate(AWRAL600_struc(n)%tgrow(AWRAL600_struc(n)%nhru))
        do i = 1, AWRAL600_struc(n)%nhru
            call ESMF_ConfigGetAttribute(LIS_config, AWRAL600_struc(n)%tgrow(i), rc=rc)
            call LIS_verify(rc, "AWRAL600 tgrow: not defined")
        enddo
    enddo
 
    ! characteristic time scale for vegetation senescence towards equilibrium for each hru
    call ESMF_ConfigFindLabel(LIS_config, "AWRAL600 tsenc:", rc = rc)
    do n=1, LIS_rc%nnest
        allocate(AWRAL600_struc(n)%tsenc(AWRAL600_struc(n)%nhru))
        do i = 1, AWRAL600_struc(n)%nhru
            call ESMF_ConfigGetAttribute(LIS_config, AWRAL600_struc(n)%tsenc(i), rc=rc)
            call LIS_verify(rc, "AWRAL600 tsenc: not defined")
        enddo
    enddo
 
    ! maximum possible root water uptake from the deep soil store for each hru
    call ESMF_ConfigFindLabel(LIS_config, "AWRAL600 ud0:", rc = rc)
    do n=1, LIS_rc%nnest
        allocate(AWRAL600_struc(n)%ud0(AWRAL600_struc(n)%nhru))
        do i = 1, AWRAL600_struc(n)%nhru
            call ESMF_ConfigGetAttribute(LIS_config, AWRAL600_struc(n)%ud0(i), rc=rc)
            call LIS_verify(rc, "AWRAL600 ud0: not defined")
        enddo
    enddo
 
    ! maximum possible root water uptake from the shallow soil store for each hru
    call ESMF_ConfigFindLabel(LIS_config, "AWRAL600 us0:", rc = rc)
    do n=1, LIS_rc%nnest
        allocate(AWRAL600_struc(n)%us0(AWRAL600_struc(n)%nhru))
        do i = 1, AWRAL600_struc(n)%nhru
            call ESMF_ConfigGetAttribute(LIS_config, AWRAL600_struc(n)%us0(i), rc=rc)
            call LIS_verify(rc, "AWRAL600 us0: not defined")
        enddo
    enddo
 
    ! vegetation photosynthetic capacity index per unit canopy cover for each hru
    call ESMF_ConfigFindLabel(LIS_config, "AWRAL600 vc:", rc = rc)
    do n=1, LIS_rc%nnest
        allocate(AWRAL600_struc(n)%vc(AWRAL600_struc(n)%nhru))
        do i = 1, AWRAL600_struc(n)%nhru
            call ESMF_ConfigGetAttribute(LIS_config, AWRAL600_struc(n)%vc(i), rc=rc)
            call LIS_verify(rc, "AWRAL600 vc: not defined")
        enddo
    enddo
 
    ! limiting the value of the relative soil moisture content of the top soil layer at which evaporation is reduced for each hru
    call ESMF_ConfigFindLabel(LIS_config, "AWRAL600 w0lime:", rc = rc)
    do n=1, LIS_rc%nnest
        allocate(AWRAL600_struc(n)%w0lime(AWRAL600_struc(n)%nhru))
        do i = 1, AWRAL600_struc(n)%nhru
            call ESMF_ConfigGetAttribute(LIS_config, AWRAL600_struc(n)%w0lime(i), rc=rc)
            call LIS_verify(rc, "AWRAL600 w0lime: not defined")
        enddo
    enddo
 
    ! Reference value of w0 that determines the rate of albedo decrease with wetness for each hru
    call ESMF_ConfigFindLabel(LIS_config, "AWRAL600 w0ref_alb:", rc = rc)
    do n=1, LIS_rc%nnest
        allocate(AWRAL600_struc(n)%w0ref_alb(AWRAL600_struc(n)%nhru))
        do i = 1, AWRAL600_struc(n)%nhru
            call ESMF_ConfigGetAttribute(LIS_config, AWRAL600_struc(n)%w0ref_alb(i), rc=rc)
            call LIS_verify(rc, "AWRAL600 w0ref_alb: not defined")
        enddo
    enddo
 
    ! water-limiting relative water content of the deep soil store for each hru
    call ESMF_ConfigFindLabel(LIS_config, "AWRAL600 wdlimu:", rc = rc)
    do n=1, LIS_rc%nnest
        allocate(AWRAL600_struc(n)%wdlimu(AWRAL600_struc(n)%nhru))
        do i = 1, AWRAL600_struc(n)%nhru
            call ESMF_ConfigGetAttribute(LIS_config, AWRAL600_struc(n)%wdlimu(i), rc=rc)
            call LIS_verify(rc, "AWRAL600 wdlimu: not defined")
        enddo
    enddo
 
    ! water-limiting relative water content of the shallow soil store for each hru
    call ESMF_ConfigFindLabel(LIS_config, "AWRAL600 wslimu:", rc = rc)
    do n=1, LIS_rc%nnest
        allocate(AWRAL600_struc(n)%wslimu(AWRAL600_struc(n)%nhru))
        do i = 1, AWRAL600_struc(n)%nhru
            call ESMF_ConfigGetAttribute(LIS_config, AWRAL600_struc(n)%wslimu(i), rc=rc)
            call LIS_verify(rc, "AWRAL600 wslimu: not defined")
        enddo
    enddo
 
    ! number of daily timesteps
    call ESMF_ConfigFindLabel(LIS_config, "AWRAL600 timesteps:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, AWRAL600_struc(n)%timesteps, rc=rc)
        call LIS_verify(rc, "AWRAL600 timesteps: not defined")
    enddo

    ! The following lines hard code the LDT NetCDF variable names. 
    do n=1, LIS_rc%nnest
        AWRAL600_struc(n)%LDT_ncvar_k_rout   = 'AWRAL600_K_ROUT'
        AWRAL600_struc(n)%LDT_ncvar_kssat    = 'AWRAL600_KSSAT'
        AWRAL600_struc(n)%LDT_ncvar_prefr    = 'AWRAL600_PREFR'
        AWRAL600_struc(n)%LDT_ncvar_s0max    = 'AWRAL600_S0MAX'
        AWRAL600_struc(n)%LDT_ncvar_slope    = 'AWRAL600_SLOPE'
        AWRAL600_struc(n)%LDT_ncvar_ssmax    = 'AWRAL600_SSMAX'
        AWRAL600_struc(n)%LDT_ncvar_k_gw     = 'AWRAL600_K_GW'
        AWRAL600_struc(n)%LDT_ncvar_kr_sd    = 'AWRAL600_KR_SD'
        AWRAL600_struc(n)%LDT_ncvar_kr_0s    = 'AWRAL600_KR_0S'
        AWRAL600_struc(n)%LDT_ncvar_k0sat    = 'AWRAL600_K0SAT'
        AWRAL600_struc(n)%LDT_ncvar_sdmax    = 'AWRAL600_SDMAX'
        AWRAL600_struc(n)%LDT_ncvar_kdsat    = 'AWRAL600_KDSAT'
        AWRAL600_struc(n)%LDT_ncvar_ne       = 'AWRAL600_NE'
        AWRAL600_struc(n)%LDT_ncvar_height   = 'AWRAL600_HEIGHT'
        AWRAL600_struc(n)%LDT_ncvar_fhru     = 'AWRAL600_FHRU'
        AWRAL600_struc(n)%LDT_ncvar_hveg     = 'AWRAL600_HVEG'
        AWRAL600_struc(n)%LDT_ncvar_laimax   = 'AWRAL600_LAIMAX'
    enddo

    ! set default restart format to netcdf
    do n=1,LIS_rc%nnest
        AWRAL600_struc(n)%rformat = "netcdf"
    enddo
    ! restart run, read restart file
    if (trim(LIS_rc%startcode) == "restart") then 
        Call ESMF_ConfigFindLabel(LIS_config, "AWRAL600 restart file:", rc=rc)
        do n=1,LIS_rc%nnest
            call ESMF_ConfigGetAttribute(LIS_config, AWRAL600_struc(n)%rfile, rc=rc)
            call LIS_verify(rc, "AWRAL600 restart file: not defined")
        enddo
        
        Call ESMF_ConfigFindLabel(LIS_config, "AWRAL600 restart file format:", rc=rc)
        do n=1,LIS_rc%nnest
            call ESMF_ConfigGetAttribute(LIS_config, AWRAL600_struc(n)%rformat, rc=rc)
            call LIS_verify(rc, "AWRAL600 restart file format: not defined")
        enddo
    ! cold start run, read initial state variables
    else 
        ! volume of water in the surface water store
        call ESMF_ConfigFindLabel(LIS_config, "AWRAL600 initial sr:", rc = rc)
        do n=1,LIS_rc%nnest
            call ESMF_ConfigGetAttribute(LIS_config, AWRAL600_struc(n)%init_sr, rc=rc)
            call LIS_verify(rc, "AWRAL600 initial sr: not defined")
        enddo

        ! groundwater storage in the unconfined aquifer
        call ESMF_ConfigFindLabel(LIS_config, "AWRAL600 initial sg:", rc = rc)
        do n=1,LIS_rc%nnest
            call ESMF_ConfigGetAttribute(LIS_config, AWRAL600_struc(n)%init_sg, rc=rc)
            call LIS_verify(rc, "AWRAL600 initial sg: not defined")
        enddo

        ! water storage in the surface soil layer for each hru
        call ESMF_ConfigFindLabel(LIS_config, "AWRAL600 initial s0:", rc = rc)
        do n=1,LIS_rc%nnest
            allocate(AWRAL600_struc(n)%init_s0(AWRAL600_struc(n)%nhru)) 
            do i=1, AWRAL600_struc(n)%nhru
                call ESMF_ConfigGetAttribute(LIS_config, AWRAL600_struc(n)%init_s0(i), rc=rc)
            end do
            call LIS_verify(rc, "AWRAL600 initial s0: not defined")
        enddo

        ! water content of the shallow soil store for each hru
        call ESMF_ConfigFindLabel(LIS_config, "AWRAL600 initial ss:", rc = rc)
        do n=1,LIS_rc%nnest
            allocate(AWRAL600_struc(n)%init_ss(AWRAL600_struc(n)%nhru))
            do i=1, AWRAL600_struc(n)%nhru
                call ESMF_ConfigGetAttribute(LIS_config, AWRAL600_struc(n)%init_ss(i), rc=rc)
            end do
            call LIS_verify(rc, "AWRAL600 initial ss: not defined")
        enddo

        ! water content of the deep soil store for each hru
        call ESMF_ConfigFindLabel(LIS_config, "AWRAL600 initial sd:", rc = rc)
        do n=1,LIS_rc%nnest
            allocate(AWRAL600_struc(n)%init_sd(AWRAL600_struc(n)%nhru))
            do i=1, AWRAL600_struc(n)%nhru
                call ESMF_ConfigGetAttribute(LIS_config, AWRAL600_struc(n)%init_sd(i), rc=rc)
            end do
            call LIS_verify(rc, "AWRAL600 initial sd: not defined")
        enddo

        ! leaf biomass
        call ESMF_ConfigFindLabel(LIS_config, "AWRAL600 initial mleaf:", rc = rc)
        do n=1,LIS_rc%nnest
            allocate(AWRAL600_struc(n)%init_mleaf(AWRAL600_struc(n)%nhru))
            do i=1, AWRAL600_struc(n)%nhru
                call ESMF_ConfigGetAttribute(LIS_config, AWRAL600_struc(n)%init_mleaf(i), rc=rc)
            end do
            call LIS_verify(rc, "AWRAL600 initial mleaf: not defined")
        enddo

    end if
     
    write(LIS_logunit, *) "[INFO] Finish reading LIS configuration file for AWRAL600 model"
     
end subroutine AWRAL600_readcrd
