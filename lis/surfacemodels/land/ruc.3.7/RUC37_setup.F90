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
! !ROUTINE: RUC37_setup
! \label{RUC37_setup}
!
! !REVISION HISTORY:
!  This subroutine is generated with the Model Implementation Toolkit developed
!  by Shugong Wang for the NASA Land Information System Version 7. The initial 
!  specification of the subroutine is defined by Sujay Kumar. 
!   1/15/15: Shugong Wang; initial implementation for LIS 7 and RUC37
!
! !INTERFACE:
subroutine RUC37_setup()
! !USES:
    use LIS_logMod,    only: LIS_logunit, LIS_verify, LIS_endrun
    use LIS_fileIOMod, only: LIS_read_param
    use LIS_coreMod,   only: LIS_rc, LIS_surface
    use RUC37_lsmMod

!
! !DESCRIPTION:
!
!  This routine is the entry point to set up the parameters
!  required for RUC37.  These include: 
!    vegetype     - vegetation category [-]
!    soiltype     - soil category [-]
!    albbck       - background snow-free albedo (0.0-1.0). [-]
!    tbot         - deep-soil time-invariant temperature (k).  representing sort of a mean annual air temperature. [K]
!    snoalb       - maximum snow albedo over deep snow (0.0-1.0) [-]
! 
!  The routines invoked are:
!  \begin{description}
!  \item[LIS\_read\_param](\ref{LIS_read_param}) \newline
!    retrieves LIS parameter data from NetCDF file
!  \item[RUC37\_read\_MULTILEVEL\_param](\ref{RUC37_read_MULTILEVEL_param}) \newline
!    retrieves MULTILEVEL spatial parameter from NetCDF file
!  \end{description}
!EOP
    implicit none
    integer           :: mtype
    integer           :: t, k, n
    integer           :: col, row
    real, allocatable :: placeholder(:,:)
    
    mtype = LIS_rc%lsm_index
    
    do n=1, LIS_rc%nnest
        ! allocate memory for place holder for #n nest
        allocate(placeholder(LIS_rc%lnc(n), LIS_rc%lnr(n)))
        
        !------------------------------------!
        ! reading spatial spatial parameters !
        !------------------------------------!
        ! read: vegetype
        if(LIS_rc%uselcmap(n) .ne. 'none') then
            write(LIS_logunit,*) "RUC37: retrieve parameter VEGETYPE from LIS"
            do t=1, LIS_rc%npatch(n, mtype)
                RUC37_struc(n)%ruc37(t)%vegetype= LIS_surface(n, mtype)%tile(t)%vegt
            enddo
        else 
            write(LIS_logunit,*) "RUC37: reading parameter VEGETYPE from ", trim(LIS_rc%paramfile(n))
            call LIS_read_param(n, trim(RUC37_struc(n)%LDT_ncvar_vegetype), placeholder)
            do t = 1, LIS_rc%npatch(n, mtype)
                col = LIS_surface(n, mtype)%tile(t)%col
                row = LIS_surface(n, mtype)%tile(t)%row
                RUC37_struc(n)%ruc37(t)%vegetype = placeholder(col, row)
            enddo 
        endif

        ! read: soiltype
        if(LIS_rc%usetexturemap(n) .ne. 'none') then
            write(LIS_logunit,*) "RUC37: retrieve parameter SOILTYPE from LIS"
            do t=1, LIS_rc%npatch(n, mtype)
                RUC37_struc(n)%ruc37(t)%soiltype= LIS_surface(n, mtype)%tile(t)%soilt
            enddo
        else 
            write(LIS_logunit,*) "RUC37: reading parameter SOILTYPE from ", trim(LIS_rc%paramfile(n))
            call LIS_read_param(n, trim(RUC37_struc(n)%LDT_ncvar_soiltype), placeholder)
            do t = 1, LIS_rc%npatch(n, mtype)
                col = LIS_surface(n, mtype)%tile(t)%col
                row = LIS_surface(n, mtype)%tile(t)%row
                RUC37_struc(n)%ruc37(t)%soiltype = placeholder(col, row)
            enddo 
        endif

        ! read: albbck
!        write(LIS_logunit,*) "RUC37: reading parameter ALBBCK from ", trim(LIS_rc%paramfile(n))
!        call LIS_read_param(n, trim(RUC37_struc(n)%LDT_ncvar_albbck), placeholder)
!        do t = 1, LIS_rc%npatch(n, mtype)
!            col = LIS_surface(n, mtype)%tile(t)%col
!            row = LIS_surface(n, mtype)%tile(t)%row
!            RUC37_struc(n)%ruc37(t)%albbck = placeholder(col, row)
!        enddo 

        ! read: tbot
        write(LIS_logunit,*) "RUC37: reading parameter TBOT from ", trim(LIS_rc%paramfile(n))
        call LIS_read_param(n, trim(RUC37_struc(n)%LDT_ncvar_tbot), placeholder)
        do t = 1, LIS_rc%npatch(n, mtype)
            col = LIS_surface(n, mtype)%tile(t)%col
            row = LIS_surface(n, mtype)%tile(t)%row
            RUC37_struc(n)%ruc37(t)%tbot = placeholder(col, row)
        enddo 

        ! read: snoalb
        write(LIS_logunit,*) "RUC37: reading parameter SNOALB from ", trim(LIS_rc%paramfile(n))
        call LIS_read_param(n, trim(RUC37_struc(n)%LDT_ncvar_snoalb), placeholder)
        do t = 1, LIS_rc%npatch(n, mtype)
            col = LIS_surface(n, mtype)%tile(t)%col
            row = LIS_surface(n, mtype)%tile(t)%row
            RUC37_struc(n)%ruc37(t)%snoalb = placeholder(col, row)
        enddo 

        !----------------------------------------------!
        ! MULTILEVEL reading spatial parameters        !
        !----------------------------------------------!
        ! read: albedo_monthly
        write(LIS_logunit,*) "RUC37: reading parameter ALBEDO_MONTHLY from ", trim(LIS_rc%paramfile(n))
        do k = 1, 12
            call RUC37_read_MULTILEVEL_param(n, RUC37_struc(n)%LDT_ncvar_albedo_monthly, k, placeholder)
            do t = 1, LIS_rc%npatch(n, mtype)
                col = LIS_surface(n, mtype)%tile(t)%col
                row = LIS_surface(n, mtype)%tile(t)%row
                RUC37_struc(n)%ruc37(t)%albedo_monthly(k) = placeholder(col, row)
            enddo 
        enddo 

        ! read: shdfac_monthly
        write(LIS_logunit,*) "RUC37: reading parameter SHDFAC_MONTHLY from ", trim(LIS_rc%paramfile(n))
        do k = 1, 12
            call RUC37_read_MULTILEVEL_param(n, RUC37_struc(n)%LDT_ncvar_shdfac_monthly, k, placeholder)
            do t = 1, LIS_rc%npatch(n, mtype)
                col = LIS_surface(n, mtype)%tile(t)%col
                row = LIS_surface(n, mtype)%tile(t)%row
                RUC37_struc(n)%ruc37(t)%shdfac_monthly(k) = placeholder(col, row)
            enddo 
        enddo 

        ! read: z0brd_monthly
!        write(LIS_logunit,*) "RUC37: reading parameter Z0BRD_MONTHLY from ", trim(LIS_rc%paramfile(n))
!        do k = 1, 12
!            call RUC37_read_MULTILEVEL_param(n, RUC37_struc(n)%LDT_ncvar_z0brd_monthly, k, placeholder)
!            do t = 1, LIS_rc%npatch(n, mtype)
!                col = LIS_surface(n, mtype)%tile(t)%col
!                row = LIS_surface(n, mtype)%tile(t)%row
!                RUC37_struc(n)%ruc37(t)%z0brd_monthly(k) = placeholder(col, row)
!            enddo 
!        enddo 

        ! read: lai_monthly
        write(LIS_logunit,*) "RUC37: reading parameter LAI_MONTHLY from ", trim(LIS_rc%paramfile(n))
        do k = 1, 12
            call RUC37_read_MULTILEVEL_param(n, RUC37_struc(n)%LDT_ncvar_lai_monthly, k, placeholder)
            do t = 1, LIS_rc%npatch(n, mtype)
                col = LIS_surface(n, mtype)%tile(t)%col
                row = LIS_surface(n, mtype)%tile(t)%row
                RUC37_struc(n)%ruc37(t)%lai_monthly(k) = placeholder(col, row)
            enddo 
        enddo 

        deallocate(placeholder)
        
        ! read soil and vegetation parameter table 
        ! subroutine ruc37_soilvegprm(mminlu, mminsl, landuse_tbl, soil_tbl, gen_tbl )
        call ruc37_soilvegprm(RUC37_struc(n)%landuse_scheme_name, &
                              RUC37_struc(n)%soil_scheme_name,    &
                              RUC37_struc(n)%landuse_tbl_name,    &
                              RUC37_struc(n)%soil_tbl_name,       &
                              RUC37_struc(n)%gen_tbl_name)
  enddo


end subroutine RUC37_setup

!BOP
!
! !ROUTINE: RUC37_read_MULTILEVEL_param
!  \label{read_MULTILEVEL_param}
!
! !REVISION HISTORY:
!  03 Sept 2004: Sujay Kumar; Initial Specification for read_laiclimo
!  30 Oct  2013: Shugong Wang; Generalization for reading MULTILEVEL spatial parameter
!
! !INTERFACE:
subroutine RUC37_read_MULTILEVEL_param(n, ncvar_name, level, placeholder)
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
    integer       :: nc, nr, t, nlevel, k
    real, pointer :: level_data(:, :, :)
    real, pointer :: level_data1(:, :, :)
    logical       :: file_exists

    inquire(file=LIS_rc%paramfile(n), exist=file_exists)
    if(file_exists) then
        write(LIS_logunit, *) 'Reading '//trim(ncvar_name)//' map for level ', level

        ! open NetCDF parameter file
        ios = nf90_open(path=trim(LIS_rc%paramfile(n)), mode=NF90_NOWRITE, ncid=nid)
        call LIS_verify(ios, 'Error in nf90_open in RUC37_read_MULTILEVEL_param')

        ! inquire the ID of east-west dimension
        ios = nf90_inq_dimid(nid, 'east_west', nc_ID)
        call LIS_verify(ios, 'Error in nf90_inq_dimid in RUC37_read_MULTILEVEL_param')

        ! inquire the ID of north-south dimension
        ios = nf90_inq_dimid(nid, 'north_south', nr_ID)
        call LIS_verify(ios, 'Error in nf90_inq_dimid in RUC37_read_MULTILEVEL_param')

        ! inquire the length of east-west dimension
        ios = nf90_inquire_dimension(nid, nc_ID, len=nc)
        call LIS_verify(ios, 'Error in nf90_inquire_dimension in RUC37_read_MULTILEVEL_param')

        ! inquire the length of north-south dimension
        ios = nf90_inquire_dimension(nid, nr_ID, len=nr)
        call LIS_verify(ios, 'Error in nf90_inquire_dimension in RUC37_read_MULTILEVEL_param')

        ! inquire the ID of parameter. 
        ios = nf90_inq_varid(nid, Trim(ncvar_name), param_ID)
        call LIS_verify(ios, trim(ncvar_name)//' field not found in the LIS param file')

        ! inquire the IDs of all dimensions. The third dimension is the level dimension
        ios = nf90_inquire_variable(nid, param_ID, dimids = dimids)
        call LIS_verify(ios, trim(ncvar_name)//' failed to inquire dimensions')

        ! inquire the length of the level dimension
        ios = nf90_inquire_dimension(nid, dimids(3), len=nlevel)
        call LIS_verify(ios, trim(ncvar_name)//' failed to inquire the length of the 3rd dimension')

        ! allocate memory
        allocate(level_data(LIS_rc%gnc(n), LIS_rc%gnr(n), nlevel))

        ! inquire the variable ID of parameter 
        ios = nf90_inq_varid(nid, trim(ncvar_name), param_ID)
        call LIS_verify(ios, trim(ncvar_name)//' field not found in the LIS param file')

        ! read parameter 
        ios = nf90_get_var(nid, param_ID, level_data)
        call LIS_verify(ios, 'Error in nf90_get_var in RUC37_read_MULTILEVEL_param')

        ! close netcdf file 
        ios = nf90_close(nid)
        call LIS_verify(ios, 'Error in nf90_close in RUC37_read_MULTILEVEL_param')

        ! grab parameter at specific level
        placeholder(:, :) = & 
             level_data(LIS_ews_halo_ind(n, LIS_localPet+1):LIS_ewe_halo_ind(n, LIS_localPet+1), &
                        LIS_nss_halo_ind(n, LIS_localPet+1):LIS_nse_halo_ind(n, LIS_localPet+1), level)

        ! free memory 
        deallocate(level_data)

    else
        write(LIS_logunit, *) 'MULTILEVEL parameter data file: ', LIS_rc%paramfile(n), ' does not exist'
        write(LIS_logunit, *) 'program stopping ...'
        call LIS_endrun
    endif
 end subroutine RUC37_read_MULTILEVEL_param
                                          

