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
! !ROUTINE: AWRAL600_setup
! \label{AWRAL600_setup}
!
! !REVISION HISTORY:
!  This subroutine is generated with the Model Implementation Toolkit developed
!  by Shugong Wang for the NASA Land Information System Version 7. The initial 
!  specification of the subroutine is defined by Sujay Kumar. 
!   12/18/18: Wendy Sharples, Shugong Wang; initial implementation for LIS 7 and AWRAL600
!
! !INTERFACE:
subroutine AWRAL600_setup()
! !USES:
    use LIS_logMod,    only: LIS_logunit, LIS_verify, LIS_endrun
    use LIS_fileIOMod, only: LIS_read_param
    use LIS_coreMod,   only: LIS_rc, LIS_surface
    use AWRAL600_lsmMod

!
! !DESCRIPTION:
!
!  This routine is the entry point to set up the parameters
!  required for AWRAL600.  These include: 
!    k_rout       - rate coefficient controlling discharge to stream [-]
!    kssat        - saturated hydraulic conductivity of shallow soil layer [mm/d]
!    prefr        - reference value for precipitation [mm]
!    s0max        - maximum storage of the surface soil layer [mm]
!    slope        - slope of the land surface [%]
!    ssmax        - maximum storage of the shallow soil layer [mm]
!    k_gw         - groundwater drainage coefficient [1/d]
!    kr_sd        - routing delay factor for the deep layer [-]
!    kr_0s        - routing delay factor for the surface layer [-]
!    k0sat        - saturated hydraulic conductivity of surface soil layer [mm/d]
!    sdmax        - maximum storage of the deep soil layer [mm]
!    kdsat        - saturated hydraulic conductivity of shallow soil layer [mm/d]
!    ne           - effective porosity [-]
!    height       - elevation of a point on the hypsometric curve [m]
!    fhru         - fraction of the cell which contains shallow and deep rooted vegetation [-]
!    hveg         - vegetation height for each hru [-]
!    laimax       - leaf area index max for each hru [-]
!    
! 
!  The routines invoked are:
!  \begin{description}
!  \item[LIS\_read\_param](\ref{LIS_read_param}) \\ 
!    retrieves LIS parameter data from NetCDF file
!  \item[AWRAL600\_read\_MULTILEVEL\_param](\ref{AWRAL600_read_MULTILEVEL_param}) \\ 
!    retrieves MULTILEVEL spatial parameter from NetCDF file
!  \end{description}
!EOP
    implicit none
    integer           :: mtype
    integer           :: t, k, n
    integer           :: col, row
    real, allocatable :: placeholder(:,:)
    real, allocatable :: placeholder2(:,:)
    
    mtype = LIS_rc%lsm_index
    
    do n=1, LIS_rc%nnest
        ! allocate memory for place holder for #n nest
        allocate(placeholder(LIS_rc%lnc(n), LIS_rc%lnr(n)))
        
        !------------------------------------!
        ! reading spatial spatial parameters !
        !------------------------------------!
        ! read: k_rout
        write(LIS_logunit,*) "[INFO] AWRAL600: reading parameter K_ROUT from ", trim(LIS_rc%paramfile(n))
        call LIS_read_param(n, trim(AWRAL600_struc(n)%LDT_ncvar_k_rout), placeholder)
        do t = 1, LIS_rc%npatch(n, mtype)
            col = LIS_surface(n, mtype)%tile(t)%col
            row = LIS_surface(n, mtype)%tile(t)%row
            AWRAL600_struc(n)%awral600(t)%k_rout = placeholder(col, row)
        enddo 

        ! read: kssat
        write(LIS_logunit,*) "[INFO] AWRAL600: reading parameter KSSAT from ", trim(LIS_rc%paramfile(n))
        call LIS_read_param(n, trim(AWRAL600_struc(n)%LDT_ncvar_kssat), placeholder)
        do t = 1, LIS_rc%npatch(n, mtype)
            col = LIS_surface(n, mtype)%tile(t)%col
            row = LIS_surface(n, mtype)%tile(t)%row
            AWRAL600_struc(n)%awral600(t)%kssat = placeholder(col, row)
        enddo 

        ! read: prefr
        write(LIS_logunit,*) "[INFO] AWRAL600: reading parameter PREFR from ", trim(LIS_rc%paramfile(n))
        call LIS_read_param(n, trim(AWRAL600_struc(n)%LDT_ncvar_prefr), placeholder)
        do t = 1, LIS_rc%npatch(n, mtype)
            col = LIS_surface(n, mtype)%tile(t)%col
            row = LIS_surface(n, mtype)%tile(t)%row
            AWRAL600_struc(n)%awral600(t)%prefr = placeholder(col, row)
        enddo 

        ! read: s0max
        write(LIS_logunit,*) "[INFO] AWRAL600: reading parameter S0MAX from ", trim(LIS_rc%paramfile(n))
        call LIS_read_param(n, trim(AWRAL600_struc(n)%LDT_ncvar_s0max), placeholder)
        do t = 1, LIS_rc%npatch(n, mtype)
            col = LIS_surface(n, mtype)%tile(t)%col
            row = LIS_surface(n, mtype)%tile(t)%row
            AWRAL600_struc(n)%awral600(t)%s0max = placeholder(col, row)
        enddo 

        ! read: slope
        write(LIS_logunit,*) "[INFO] AWRAL600: reading parameter SLOPE from ", trim(LIS_rc%paramfile(n))
        call LIS_read_param(n, trim(AWRAL600_struc(n)%LDT_ncvar_slope), placeholder)
        do t = 1, LIS_rc%npatch(n, mtype)
            col = LIS_surface(n, mtype)%tile(t)%col
            row = LIS_surface(n, mtype)%tile(t)%row
            AWRAL600_struc(n)%awral600(t)%slope = placeholder(col, row)
        enddo 

        ! read: ssmax
        write(LIS_logunit,*) "[INFO] AWRAL600: reading parameter SSMAX from ", trim(LIS_rc%paramfile(n))
        call LIS_read_param(n, trim(AWRAL600_struc(n)%LDT_ncvar_ssmax), placeholder)
        do t = 1, LIS_rc%npatch(n, mtype)
            col = LIS_surface(n, mtype)%tile(t)%col
            row = LIS_surface(n, mtype)%tile(t)%row
            AWRAL600_struc(n)%awral600(t)%ssmax = placeholder(col, row)
        enddo 

        ! read: k_gw
        write(LIS_logunit,*) "[INFO] AWRAL600: reading parameter K_GW from ", trim(LIS_rc%paramfile(n))
        call LIS_read_param(n, trim(AWRAL600_struc(n)%LDT_ncvar_k_gw), placeholder)
        do t = 1, LIS_rc%npatch(n, mtype)
            col = LIS_surface(n, mtype)%tile(t)%col
            row = LIS_surface(n, mtype)%tile(t)%row
            AWRAL600_struc(n)%awral600(t)%k_gw = placeholder(col, row)
        enddo 

        ! read: kr_sd
        write(LIS_logunit,*) "[INFO] AWRAL600: reading parameter KR_SD from ", trim(LIS_rc%paramfile(n))
        call LIS_read_param(n, trim(AWRAL600_struc(n)%LDT_ncvar_kr_sd), placeholder)
        do t = 1, LIS_rc%npatch(n, mtype)
            col = LIS_surface(n, mtype)%tile(t)%col
            row = LIS_surface(n, mtype)%tile(t)%row
            AWRAL600_struc(n)%awral600(t)%kr_sd = placeholder(col, row)
        enddo 

        ! read: kr_0s
        write(LIS_logunit,*) "[INFO] AWRAL600: reading parameter KR_0S from ", trim(LIS_rc%paramfile(n))
        call LIS_read_param(n, trim(AWRAL600_struc(n)%LDT_ncvar_kr_0s), placeholder)
        do t = 1, LIS_rc%npatch(n, mtype)
            col = LIS_surface(n, mtype)%tile(t)%col
            row = LIS_surface(n, mtype)%tile(t)%row
            AWRAL600_struc(n)%awral600(t)%kr_0s = placeholder(col, row)
        enddo 

        ! read: k0sat
        write(LIS_logunit,*) "[INFO] AWRAL600: reading parameter K0SAT from ", trim(LIS_rc%paramfile(n))
        call LIS_read_param(n, trim(AWRAL600_struc(n)%LDT_ncvar_k0sat), placeholder)
        do t = 1, LIS_rc%npatch(n, mtype)
            col = LIS_surface(n, mtype)%tile(t)%col
            row = LIS_surface(n, mtype)%tile(t)%row
            AWRAL600_struc(n)%awral600(t)%k0sat = placeholder(col, row)
        enddo 

        ! read: sdmax
        write(LIS_logunit,*) "[INFO] AWRAL600: reading parameter SDMAX from ", trim(LIS_rc%paramfile(n))
        call LIS_read_param(n, trim(AWRAL600_struc(n)%LDT_ncvar_sdmax), placeholder)
        do t = 1, LIS_rc%npatch(n, mtype)
            col = LIS_surface(n, mtype)%tile(t)%col
            row = LIS_surface(n, mtype)%tile(t)%row
            AWRAL600_struc(n)%awral600(t)%sdmax = placeholder(col, row)
        enddo 

        ! read: kdsat
        write(LIS_logunit,*) "[INFO] AWRAL600: reading parameter KDSAT from ", trim(LIS_rc%paramfile(n))
        call LIS_read_param(n, trim(AWRAL600_struc(n)%LDT_ncvar_kdsat), placeholder)
        do t = 1, LIS_rc%npatch(n, mtype)
            col = LIS_surface(n, mtype)%tile(t)%col
            row = LIS_surface(n, mtype)%tile(t)%row
            AWRAL600_struc(n)%awral600(t)%kdsat = placeholder(col, row)
        enddo 

        ! read: ne
        write(LIS_logunit,*) "[INFO] AWRAL600: reading parameter NE from ", trim(LIS_rc%paramfile(n))
        call LIS_read_param(n, trim(AWRAL600_struc(n)%LDT_ncvar_ne), placeholder)
        do t = 1, LIS_rc%npatch(n, mtype)
            col = LIS_surface(n, mtype)%tile(t)%col
            row = LIS_surface(n, mtype)%tile(t)%row
            AWRAL600_struc(n)%awral600(t)%ne = placeholder(col, row)
        enddo 

        !----------------------------------------------!
        ! MULTILEVEL reading spatial spatial parameters !
        !----------------------------------------------!
        ! read: height
        write(LIS_logunit,*) "[INFO] AWRAL600: reading parameter HEIGHT from ", trim(LIS_rc%paramfile(n))
        do k = 1, AWRAL600_struc(n)%nhypsbins
            call AWRAL600_read_MULTILEVEL_param(n, AWRAL600_struc(n)%LDT_ncvar_height, k, placeholder)
            do t = 1, LIS_rc%npatch(n, mtype)
                col = LIS_surface(n, mtype)%tile(t)%col
                row = LIS_surface(n, mtype)%tile(t)%row
                AWRAL600_struc(n)%awral600(t)%height(k) = placeholder(col, row)
            enddo 
        enddo 

        ! read: fhru
        write(LIS_logunit,*) "[INFO] AWRAL600: reading parameter FHRU from ", trim(LIS_rc%paramfile(n))
        do k = 1, AWRAL600_struc(n)%nhru
            call AWRAL600_read_MULTILEVEL_param(n, AWRAL600_struc(n)%LDT_ncvar_fhru, k, placeholder)
            do t = 1, LIS_rc%npatch(n, mtype)
                col = LIS_surface(n, mtype)%tile(t)%col
                row = LIS_surface(n, mtype)%tile(t)%row
                AWRAL600_struc(n)%awral600(t)%fhru(k) = placeholder(col, row)
            enddo 
        enddo 

        ! read: hveg
        write(LIS_logunit,*) "[INFO] AWRAL600: reading parameter HVEG from ", trim(LIS_rc%paramfile(n))
        do k = 1, AWRAL600_struc(n)%nhru
            call AWRAL600_read_MULTILEVEL_param(n, AWRAL600_struc(n)%LDT_ncvar_hveg, k, placeholder)
            do t = 1, LIS_rc%npatch(n, mtype)
                col = LIS_surface(n, mtype)%tile(t)%col
                row = LIS_surface(n, mtype)%tile(t)%row
                AWRAL600_struc(n)%awral600(t)%hveg(k) = placeholder(col, row)
            enddo 
        enddo 

        ! read: laimax
        write(LIS_logunit,*) "[INFO] AWRAL600: reading parameter LAIMAX from ", trim(LIS_rc%paramfile(n))
        do k = 1, AWRAL600_struc(n)%nhru
            call AWRAL600_read_MULTILEVEL_param(n, AWRAL600_struc(n)%LDT_ncvar_laimax, k, placeholder)
            do t = 1, LIS_rc%npatch(n, mtype)
                col = LIS_surface(n, mtype)%tile(t)%col
                row = LIS_surface(n, mtype)%tile(t)%row
                AWRAL600_struc(n)%awral600(t)%laimax(k) = placeholder(col, row)
            enddo 
        enddo
      
        deallocate(placeholder)
    enddo


end subroutine AWRAL600_setup

!BOP
!
! !ROUTINE: AWRAL600_read_MULTILEVEL_param
!  \label{read_MULTILEVEL_param}
!
! !INTERFACE:
subroutine AWRAL600_read_MULTILEVEL_param(n, ncvar_name, level, placeholder)
! !USES:
    use netcdf
    use LIS_coreMod, only : LIS_rc, LIS_domain, LIS_localPet,   &   
                            LIS_ews_halo_ind, LIS_ewe_halo_ind, &
                            LIS_nss_halo_ind, LIS_nse_halo_ind   
    use LIS_logMod,  only : LIS_logunit, LIS_verify, LIS_endrun
    use LIS_fileIOMod, only: LIS_read_param
    implicit none
! !ARGUMENTS: 
    integer, intent(in)          :: n
    integer, intent(in)          :: level
    character(len=*), intent(in) :: ncvar_name 
    real, intent(out)            :: placeholder(LIS_rc%lnc(n), LIS_rc%lnr(n))
! !DESCRIPTION:
!  This subroutine reads MULTILEVEL parameters from the LIS
!  NetCDF parameter data file
!  
!  The arguments are:
!  \begin{description}
!   \item[n]
!    index of n
!   \item[level]
!    level index (month, quarter, soil layer, snow layer) of the data to be read
!   \item[array]
!    array containing returned values
!   \end{description}
!
!EOP      

    integer       :: ios1
    integer       :: ios, nid, param_ID, nc_ID, nr_ID, dimids(3)
    integer       :: nc, nr, t, k
    integer       :: nlevel
    real, allocatable :: level_data(:, :, :)
    logical       :: file_exists

    inquire(file=LIS_rc%paramfile(n), exist=file_exists)
    if(file_exists) then
        write(LIS_logunit, *) '[INFO] Reading '//trim(ncvar_name)//' map for level ', level

        ! open NetCDF parameter file
        ios = nf90_open(path=trim(LIS_rc%paramfile(n)), mode=NF90_NOWRITE, ncid=nid)
        call LIS_verify(ios, 'Error in nf90_open in AWRAL600_read_MULTILEVEL_param')

        ! inquire the ID of east-west dimension
        ios = nf90_inq_dimid(nid, 'east_west', nc_ID)
        call LIS_verify(ios, 'Error in nf90_inq_dimid in AWRAL600_read_MULTILEVEL_param')

        ! inquire the ID of north-south dimension
        ios = nf90_inq_dimid(nid, 'north_south', nr_ID)
        call LIS_verify(ios, 'Error in nf90_inq_dimid in AWRAL600_read_MULTILEVEL_param')

        ! inquire the length of east-west dimension
        ios = nf90_inquire_dimension(nid, nc_ID, len=nc)
        call LIS_verify(ios, 'Error in nf90_inquire_dimension in AWRAL600_read_MULTILEVEL_param')

        ! inquire the length of north-south dimension
        ios = nf90_inquire_dimension(nid, nr_ID, len=nr)
        call LIS_verify(ios, 'Error in nf90_inquire_dimension in AWRAL600_read_MULTILEVEL_param')

        ! inquire the ID of parameter. 
        ios = nf90_inq_varid(nid, Trim(ncvar_name), param_ID)
        call LIS_verify(ios, trim(ncvar_name)//'field not found in the LIS param file')

        ! inquire the IDs of all dimensions. The third dimension is the level dimension
        ios = nf90_inquire_variable(nid, param_ID, dimids = dimids)
        call LIS_verify(ios, trim(ncvar_name)//'failed to inquire dimensions')

        ! inquire the length of the level dimension
        ios = nf90_inquire_dimension(nid, dimids(3), len=nlevel)
        call LIS_verify(ios, trim(ncvar_name)//'failed to inquire the length of the 3rd dimension')

        ! allocate memory
        allocate(level_data (LIS_rc%gnc(n), LIS_rc%gnr(n), nlevel))

        ! inquire the variable ID of parameter 
        ios = nf90_inq_varid(nid, trim(ncvar_name), param_ID)
        call LIS_verify(ios, trim(ncvar_name)//'field not found in the LIS param file')

        ! read parameter 
        ios = nf90_get_var(nid, param_ID, level_data)
        call LIS_verify(ios, 'Error in nf90_get_var in AWRAL600_read_MULTILEVEL_param')

        ! close netcdf file 
        ios = nf90_close(nid)
        call LIS_verify(ios, 'Error in nf90_close in AWRAL600_read_MULTILEVEL_param')

        ! grab parameter at specific level
        placeholder(:, :) = & 
             level_data(LIS_ews_halo_ind(n, LIS_localPet+1):LIS_ewe_halo_ind(n, LIS_localPet+1), &
                        LIS_nss_halo_ind(n, LIS_localPet+1):LIS_nse_halo_ind(n, LIS_localPet+1), level)

        ! free memory 
        deallocate(level_data)

    else
        write(LIS_logunit, *) '[ERR] MULTILEVEL parameter data file: ', LIS_rc%paramfile(n), ' does not exist'
        write(LIS_logunit, *) '[ERR] program stopping ...'
        call LIS_endrun
    endif
 end subroutine AWRAL600_read_MULTILEVEL_param
                                          

