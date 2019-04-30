!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.0     
!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
#include "LIS_misc.h"
!BOP
!
! !ROUTINE: NoahMP401_setup
! \label{NoahMP401_setup}
!
! !REVISION HISTORY:
!  This subroutine is generated with the Model Implementation Toolkit developed
!  by Shugong Wang for the NASA Land Information System Version 7. The initial 
!  specification of the subroutine is defined by Sujay Kumar. 
!   10/25/18: Shugong Wang, Zhuo Wang; initial implementation for LIS 7 and NoahMP401
!
! !INTERFACE:
subroutine NoahMP401_setup()
! !USES:
    use NoahMP401_lsmMod
    use LIS_logMod,    only: LIS_logunit, LIS_verify, LIS_endrun
    use LIS_fileIOMod, only: LIS_read_param!, LIS_convertParamDataToLocalDomain
    use LIS_coreMod,   only: LIS_rc, LIS_surface
    use NOAHMP_TABLES_401, only: read_mp_veg_parameters,     &
                                 read_mp_soil_parameters,    &
                                 read_mp_rad_parameters,     &
                                 read_mp_global_parameters,  &
                                 read_mp_crop_parameters,    &
                                 read_mp_optional_parameters
   
!
! !DESCRIPTION:
!
!  This routine is the entry point to set up the parameters
!  required for NoahMP401.  These include: 
!    vegetype     - vegetation type [-]
!    soiltype     - soil type [-]
!    tbot         - deep soil temperature [K]
!    planting     - planting date [-]
!    harvest      - harvest date [-]
!    season_gdd   - growing season GDD [-]
!    soilcL1      - soil texture in layer 1 [-]
!    soilcL2      - soil texture in layer 2 [-]
!    soilcL3      - soil texture in layer 3 [-]
!    soilcL4      - soil texture in layer 4 [-]
! 
!  The routines invoked are:
!  \begin{description}
!  \item[LIS\_read\_param](\ref{LIS_read_param}) \\ 
!    retrieves LIS parameter data from NetCDF file
!  \item[NOAHMP401\_read\_MULTILEVEL\_param](\ref{NOAHMP401_read_MULTILEVEL_param}) \\ 
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
        
        !----------------------------!
        ! reading spatial parameters !
        !----------------------------!
        ! vegetype takes value from the LIS built-in parameter vegt
        !TODO: convert vegetation data source into vegetation types
        if(LIS_rc%uselcmap(n) .ne. 'none') then
            write(LIS_logunit,*) &
             "[INFO] Noah-MP.4.0.1 retrieve parameter VEGETYPE from LIS"
            do t=1, LIS_rc%npatch(n, mtype)
                NOAHMP401_struc(n)%noahmp401(t)%vegetype= LIS_surface(n, mtype)%tile(t)%vegt
            enddo
        else 
            ! read: vegetype
            write(LIS_logunit,*) &
             "[INFO] Noah-MP.4.0.1 reading parameter VEGETYPE from ", &
                                   trim(LIS_rc%paramfile(n))
            call LIS_read_param(n, trim(NOAHMP401_struc(n)%LDT_ncvar_vegetype), placeholder)
            do t = 1, LIS_rc%npatch(n, mtype)
                col = LIS_surface(n, mtype)%tile(t)%col
                row = LIS_surface(n, mtype)%tile(t)%row
                NOAHMP401_struc(n)%noahmp401(t)%vegetype = placeholder(col, row)
            enddo
        endif
        ! soiltype takes value from the LIS built-in parameter soilt
        !TODO: convert soil texture into soil types according to scheme
        if(LIS_rc%usetexturemap(n) .ne. 'none') then
            write(LIS_logunit,*) &
             "[INFO] Noah-MP.4.0.1 retrieve parameter SOILTYPE from LIS"
            do t=1, LIS_rc%npatch(n, mtype)
                NOAHMP401_struc(n)%noahmp401(t)%soiltype= LIS_surface(n, mtype)%tile(t)%soilt
            enddo
        else 
            ! read: soiltype
            write(LIS_logunit,*) &
             "[INFO] Noah-MP.4.0.1 reading parameter SOILTYPE from ", &
                                   trim(LIS_rc%paramfile(n))
            call LIS_read_param(n, trim(NOAHMP401_struc(n)%LDT_ncvar_soiltype), placeholder)
            do t = 1, LIS_rc%npatch(n, mtype)
                col = LIS_surface(n, mtype)%tile(t)%col
                row = LIS_surface(n, mtype)%tile(t)%row
                NOAHMP401_struc(n)%noahmp401(t)%soiltype = placeholder(col, row)
            enddo 
        endif
        ! read: tbot
        write(LIS_logunit,*) &
         "[INFO] Noah-MP.4.0.1 reading parameter TBOT from ", &
                               trim(LIS_rc%paramfile(n))
        call LIS_read_param(n, trim(NOAHMP401_struc(n)%LDT_ncvar_tbot), placeholder)
        do t = 1, LIS_rc%npatch(n, mtype)
            col = LIS_surface(n, mtype)%tile(t)%col
            row = LIS_surface(n, mtype)%tile(t)%row
            NOAHMP401_struc(n)%noahmp401(t)%tbot = placeholder(col, row)
        enddo

        !!! SW 11/06/2018 
        if(NOAHMP401_struc(n)%crop_opt .ne.0) then 
            ! read: planting
            write(LIS_logunit,*) &
             "[INFO] Noah-MP.4.0.1 reading parameter PLANTING from ", trim(LIS_rc%paramfile(n))
            call LIS_read_param(n, trim(NOAHMP401_struc(n)%LDT_ncvar_planting), placeholder)
            do t = 1, LIS_rc%npatch(n, mtype)
                col = LIS_surface(n, mtype)%tile(t)%col
                row = LIS_surface(n, mtype)%tile(t)%row
                NOAHMP401_struc(n)%noahmp401(t)%planting = placeholder(col, row)
            enddo 

            ! read: harvest
            write(LIS_logunit,*) "[INFO] Noah-MP.4.0.1 reading parameter HARVEST from ", trim(LIS_rc%paramfile(n))
            call LIS_read_param(n, trim(NOAHMP401_struc(n)%LDT_ncvar_harvest), placeholder)
            do t = 1, LIS_rc%npatch(n, mtype)
                col = LIS_surface(n, mtype)%tile(t)%col
                row = LIS_surface(n, mtype)%tile(t)%row
                NOAHMP401_struc(n)%noahmp401(t)%harvest = placeholder(col, row)
            enddo 

            ! read: season_gdd
            write(LIS_logunit,*) "[INFO] Noah-MP.4.0.1 reading parameter SEASON_GDD from ", trim(LIS_rc%paramfile(n))
            call LIS_read_param(n, trim(NOAHMP401_struc(n)%LDT_ncvar_season_gdd), placeholder)
            do t = 1, LIS_rc%npatch(n, mtype)
                col = LIS_surface(n, mtype)%tile(t)%col
                row = LIS_surface(n, mtype)%tile(t)%row
                NOAHMP401_struc(n)%noahmp401(t)%season_gdd = placeholder(col, row)
            enddo 
        endif
        
        !!! SW 11/06/2018 
        if(NOAHMP401_struc(n)%soil_opt .eq. 2) then 
            ! read: soilcL1
            write(LIS_logunit,*) "[INFO] Noah-MP.4.0.1 reading parameter SOILCL1 from ", trim(LIS_rc%paramfile(n))
            call LIS_read_param(n, trim(NOAHMP401_struc(n)%LDT_ncvar_soilcL1), placeholder)
            do t = 1, LIS_rc%npatch(n, mtype)
                col = LIS_surface(n, mtype)%tile(t)%col
                row = LIS_surface(n, mtype)%tile(t)%row
                NOAHMP401_struc(n)%noahmp401(t)%soilcl1 = placeholder(col, row)
            enddo 

            ! read: soilcL2
            write(LIS_logunit,*) "[INFO] Noah-MP.4.0.1 reading parameter SOILCL2 from ", trim(LIS_rc%paramfile(n))
            call LIS_read_param(n, trim(NOAHMP401_struc(n)%LDT_ncvar_soilcL2), placeholder)
            do t = 1, LIS_rc%npatch(n, mtype)
                col = LIS_surface(n, mtype)%tile(t)%col
                row = LIS_surface(n, mtype)%tile(t)%row
                NOAHMP401_struc(n)%noahmp401(t)%soilcl2 = placeholder(col, row)
            enddo 

            ! read: soilcL3
            write(LIS_logunit,*) "[INFO] Noah-MP.4.0.1 reading parameter SOILCL3 from ", trim(LIS_rc%paramfile(n))
            call LIS_read_param(n, trim(NOAHMP401_struc(n)%LDT_ncvar_soilcL3), placeholder)
            do t = 1, LIS_rc%npatch(n, mtype)
                col = LIS_surface(n, mtype)%tile(t)%col
                row = LIS_surface(n, mtype)%tile(t)%row
                NOAHMP401_struc(n)%noahmp401(t)%soilcl3 = placeholder(col, row)
            enddo 

            ! read: soilcL4
            write(LIS_logunit,*) "[INFO] Noah-MP.4.0.1 reading parameter SOILCL4 from ", trim(LIS_rc%paramfile(n))
            call LIS_read_param(n, trim(NOAHMP401_struc(n)%LDT_ncvar_soilcL4), placeholder)
            do t = 1, LIS_rc%npatch(n, mtype)
                col = LIS_surface(n, mtype)%tile(t)%col
                row = LIS_surface(n, mtype)%tile(t)%row
                NOAHMP401_struc(n)%noahmp401(t)%soilcl4 = placeholder(col, row)
            enddo 
        endif 

        !----------------------------------------------!
        ! MULTILEVEL reading spatial spatial parameters !
        !----------------------------------------------!
        write(LIS_logunit,*) "[INFO] Noah-MP.4.0.1 reading parameter SHDFAC_MONTHLY from ", trim(LIS_rc%paramfile(n))
        do k = 1, 12
            call NOAHMP401_read_MULTILEVEL_param(n, NOAHMP401_struc(n)%LDT_ncvar_shdfac_monthly, k, placeholder)
            do t = 1, LIS_rc%npatch(n, mtype)
                col = LIS_surface(n, mtype)%tile(t)%col
                row = LIS_surface(n, mtype)%tile(t)%row
                NOAHMP401_struc(n)%noahmp401(t)%shdfac_monthly(k) = placeholder(col, row)
            enddo 
        enddo 
        
        if(NOAHMP401_struc(n)%soil_opt .eq. 3) then
            ! read: soilcomp
            write(LIS_logunit,*) "[INFO] Noah-MP.4.0.1 reading parameter SOILCOMP from ", trim(LIS_rc%paramfile(n))
            do k = 1, 8
                call NOAHMP401_read_MULTILEVEL_param(n, NOAHMP401_struc(n)%LDT_ncvar_soilcomp, k, placeholder)
                do t = 1, LIS_rc%npatch(n, mtype)
                    col = LIS_surface(n, mtype)%tile(t)%col
                    row = LIS_surface(n, mtype)%tile(t)%row
                    NOAHMP401_struc(n)%noahmp401(t)%soilcomp(k) = placeholder(col, row)
                enddo 
            enddo 
        endif
        deallocate(placeholder)

    !!!! read Noah-MP parameter tables  Shugong Wang 11/06/2018 
        write(LIS_logunit,*) "[INFO] Noah-MP.4.0.1 soil parameter table: ",     trim(NOAHMP401_struc(n)%soil_tbl_name)       
        write(LIS_logunit,*) "[INFO] Noah-MP.4.0.1 general parameter table: ",  trim(NOAHMP401_struc(n)%gen_tbl_name)        
        write(LIS_logunit,*) "[INFO] Noah-MP.4.0.1 MP parameter table: ",              trim(NOAHMP401_struc(n)%noahmp_tbl_name)     
        write(LIS_logunit,*) "[INFO] Noah-MP.4.0.1 Landuse classification scheme: ",       trim(NOAHMP401_struc(n)%landuse_scheme_name) 
        write(LIS_logunit,*) "[INFO] Noah-MP.4.0.1 Soil classification scheme: ",          "STAS (default, cannot change)" 
        call read_mp_veg_parameters(trim(NOAHMP401_struc(n)%landuse_scheme_name), trim(NOAHMP401_struc(n)%noahmp_tbl_name))
        call read_mp_soil_parameters(trim(NOAHMP401_struc(n)%soil_tbl_name), trim(NOAHMP401_struc(n)%gen_tbl_name))
        call read_mp_rad_parameters(trim(NOAHMP401_struc(n)%noahmp_tbl_name))
        call read_mp_global_parameters(trim(NOAHMP401_struc(n)%noahmp_tbl_name))
        call read_mp_crop_parameters(trim(NOAHMP401_struc(n)%noahmp_tbl_name))
        call read_mp_optional_parameters(trim(NOAHMP401_struc(n)%noahmp_tbl_name))
    enddo

end subroutine NoahMP401_setup

!BOP
!
! !ROUTINE: NOAHMP401_read_MULTILEVEL_param
!  \label{read_MULTILEVEL_param}
!
! !REVISION HISTORY:
!  03 Sept 2004: Sujay Kumar; Initial Specification for read_laiclimo
!  30 Oct  2013: Shugong Wang; Generalization for reading MULTILEVEL spatial parameter
!
! !INTERFACE:
subroutine NOAHMP401_read_MULTILEVEL_param(n, ncvar_name, level, placeholder)
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
    logical       :: file_exists

    inquire(file=LIS_rc%paramfile(n), exist=file_exists)
    if(file_exists) then
        write(LIS_logunit, *) '[INFO] Reading '//trim(ncvar_name)//' map for level ', level

        ! open NetCDF parameter file
        ios = nf90_open(path=trim(LIS_rc%paramfile(n)), mode=NF90_NOWRITE, ncid=nid)
        call LIS_verify(ios, 'Error in nf90_open in NOAHMP401_read_MULTILEVEL_param')

        ! inquire the ID of east-west dimension
        ios = nf90_inq_dimid(nid, 'east_west', nc_ID)
        call LIS_verify(ios, 'Error in nf90_inq_dimid in NOAHMP401_read_MULTILEVEL_param')

        ! inquire the ID of north-south dimension
        ios = nf90_inq_dimid(nid, 'north_south', nr_ID)
        call LIS_verify(ios, 'Error in nf90_inq_dimid in NOAHMP401_read_MULTILEVEL_param')

        ! inquire the length of east-west dimension
        ios = nf90_inquire_dimension(nid, nc_ID, len=nc)
        call LIS_verify(ios, 'Error in nf90_inquire_dimension in NOAHMP401_read_MULTILEVEL_param')

        ! inquire the length of north-south dimension
        ios = nf90_inquire_dimension(nid, nr_ID, len=nr)
        call LIS_verify(ios, 'Error in nf90_inquire_dimension in NOAHMP401_read_MULTILEVEL_param')

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
        allocate(level_data (LIS_rc%gnc(n), LIS_rc%gnr(n), nlevel))

        ! inquire the variable ID of parameter 
        ios = nf90_inq_varid(nid, trim(ncvar_name), param_ID)
        call LIS_verify(ios, trim(ncvar_name)//' field not found in the LIS param file')

        ! read parameter 
        ios = nf90_get_var(nid, param_ID, level_data)
        call LIS_verify(ios, 'Error in nf90_get_var in NOAHMP401_read_MULTILEVEL_param')

        ! close netcdf file 
        ios = nf90_close(nid)
        call LIS_verify(ios, 'Error in nf90_close in NOAHMP401_read_MULTILEVEL_param')

        ! grab parameter at specific level
        placeholder(:, :) = & 
             level_data(LIS_ews_halo_ind(n, LIS_localPet+1):LIS_ewe_halo_ind(n, LIS_localPet+1), &
                        LIS_nss_halo_ind(n, LIS_localPet+1):LIS_nse_halo_ind(n, LIS_localPet+1), level)

        ! free memory 
        deallocate(level_data)

    else
        write(LIS_logunit, *) '[ERR] MULTILEVEL parameter data file: ', &
                                     trim(LIS_rc%paramfile(n))
        write(LIS_logunit, *) '[ERR] does not exist.'
        write(LIS_logunit, *) '[ERR] program stopping ...'
        call LIS_endrun
    endif
 end subroutine NOAHMP401_read_MULTILEVEL_param
                                          

