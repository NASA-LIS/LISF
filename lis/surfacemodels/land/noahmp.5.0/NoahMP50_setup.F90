!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.4
!
! Copyright (c) 2022 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------

#include "LIS_misc.h"
!BOP
!
! !ROUTINE: NoahMP50_setup
! \label{NoahMP50_setup}
!
! !REVISION HISTORY:
!  05/01/23: Cenlin He; update to work with refactored Noah-MP (v5.0 and later)
!
! !INTERFACE:

subroutine NoahMP50_setup()

! !USES:
    use NoahMP50_lsmMod
    use NoahmpIOVarType
    use LIS_logMod,    only: LIS_logunit, LIS_verify, LIS_endrun
    use LIS_fileIOMod, only: LIS_read_param!, LIS_convertParamDataToLocalDomain
    use LIS_coreMod,   only: LIS_rc, LIS_surface
    use NoahmpReadTableMod, only : NoahmpReadTable
    use NoahmpIOVarInitMod, only : NoahmpIOVarInitDefault

!
! !DESCRIPTION:
!
!  This routine is the entry point to set up the parameters
!  required for NoahMP50.  These include: 
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
!    irfract      - total irrigation fraction [-]
!    sifract      - sprinker irrigation fraction [-]
!    mifract      - micro irrigation fraction [-]
!    fifract      - flood irrigation fraction [-]
!    tdfract      - tile drainage fraction [-]
!  The routines invoked are:
!  \begin{description}
!  \item[LIS\_read\_param](\ref{LIS_read_param}) \\ 
!    retrieves LIS parameter data from NetCDF file
!  \item[NoahMP50\_read\_MULTILEVEL\_param](\ref{NoahMP50_read_MULTILEVEL_param}) \\ 
!    retrieves MULTILEVEL spatial parameter from NetCDF file
!  \end{description}
!EOP

    implicit none
    integer           :: mtype
    integer           :: t, k, n
    integer           :: col, row
    real, allocatable :: placeholder(:,:)
    integer           :: soilcolor, vegtyp, soiltyp(4), slopetyp, croptype    
    mtype = LIS_rc%lsm_index
    
    do n=1, LIS_rc%nnest
        ! allocate memory for place holder for #n nest
        allocate(placeholder(LIS_rc%lnc(n), LIS_rc%lnr(n)))
        
        !----------------------------!
        ! reading spatial parameters !
        !----------------------------!
        ! vegetype takes value from the LIS built-in parameter vegt
        if(LIS_rc%uselcmap(n) .ne. 'none') then
            write(LIS_logunit,*) &
             "[INFO] Noah-MP.5.0 retrieve parameter VEGETYPE from LIS"
            do t=1, LIS_rc%npatch(n, mtype)
                NoahMP50_struc(n)%noahmp50(t)%vegetype= LIS_surface(n, mtype)%tile(t)%vegt
            enddo
        else 
            ! read: vegetype
            write(LIS_logunit,*) &
             "[INFO] Noah-MP.5.0 reading parameter VEGETYPE from ", &
                                   trim(LIS_rc%paramfile(n))
            call LIS_read_param(n, trim(NoahMP50_struc(n)%LDT_ncvar_vegetype), placeholder)
            do t = 1, LIS_rc%npatch(n, mtype)
                col = LIS_surface(n, mtype)%tile(t)%col
                row = LIS_surface(n, mtype)%tile(t)%row
                NoahMP50_struc(n)%noahmp50(t)%vegetype = placeholder(col, row)
            enddo
        endif
        ! soiltype takes value from the LIS built-in parameter soilt
        if(LIS_rc%usetexturemap(n) .ne. 'none') then
            write(LIS_logunit,*) &
             "[INFO] Noah-MP.5.0 retrieve parameter SOILTYPE from LIS"
            do t=1, LIS_rc%npatch(n, mtype)
                NoahMP50_struc(n)%noahmp50(t)%soiltype= LIS_surface(n, mtype)%tile(t)%soilt
            enddo
        else 
            ! read: soiltype
            write(LIS_logunit,*) &
             "[INFO] Noah-MP.5.0 reading parameter SOILTYPE from ", &
                                   trim(LIS_rc%paramfile(n))
            call LIS_read_param(n, trim(NoahMP50_struc(n)%LDT_ncvar_soiltype), placeholder)
            do t = 1, LIS_rc%npatch(n, mtype)
                col = LIS_surface(n, mtype)%tile(t)%col
                row = LIS_surface(n, mtype)%tile(t)%row
                NoahMP50_struc(n)%noahmp50(t)%soiltype = placeholder(col, row)
            enddo 
        endif
        ! read: tbot
        write(LIS_logunit,*) &
         "[INFO] Noah-MP.5.0 reading parameter TBOT from ", &
                               trim(LIS_rc%paramfile(n))
        call LIS_read_param(n, trim(NoahMP50_struc(n)%LDT_ncvar_tbot), placeholder)
        do t = 1, LIS_rc%npatch(n, mtype)
            col = LIS_surface(n, mtype)%tile(t)%col
            row = LIS_surface(n, mtype)%tile(t)%row
            NoahMP50_struc(n)%noahmp50(t)%tbot = placeholder(col, row)
        enddo

        if(NoahMP50_struc(n)%crop_opt > 0) then 
            ! read: planting
            write(LIS_logunit,*) &
             "[INFO] Noah-MP.5.0 reading parameter PLANTING from ", &
             trim(LIS_rc%paramfile(n))
            call LIS_read_param(n, trim(NoahMP50_struc(n)%LDT_ncvar_planting), placeholder)
            do t = 1, LIS_rc%npatch(n, mtype)
                col = LIS_surface(n, mtype)%tile(t)%col
                row = LIS_surface(n, mtype)%tile(t)%row
                NoahMP50_struc(n)%noahmp50(t)%planting = placeholder(col, row)
            enddo 

            ! read: harvest
            write(LIS_logunit,*) "[INFO] Noah-MP.5.0 reading parameter HARVEST from ",&
                 trim(LIS_rc%paramfile(n))
            call LIS_read_param(n, trim(NoahMP50_struc(n)%LDT_ncvar_harvest), placeholder)
            do t = 1, LIS_rc%npatch(n, mtype)
                col = LIS_surface(n, mtype)%tile(t)%col
                row = LIS_surface(n, mtype)%tile(t)%row
                NoahMP50_struc(n)%noahmp50(t)%harvest = placeholder(col, row)
            enddo 

            ! read: season_gdd
            write(LIS_logunit,*) "[INFO] Noah-MP.5.0 reading parameter SEASON_GDD from ", &
                 trim(LIS_rc%paramfile(n))
            call LIS_read_param(n, trim(NoahMP50_struc(n)%LDT_ncvar_season_gdd), placeholder)
            do t = 1, LIS_rc%npatch(n, mtype)
                col = LIS_surface(n, mtype)%tile(t)%col
                row = LIS_surface(n, mtype)%tile(t)%row
                NoahMP50_struc(n)%noahmp50(t)%season_gdd = placeholder(col, row)
            enddo 
        endif

        ! for irrigation
        if(NoahMP50_struc(n)%irr_opt > 0) then
            ! read: total irrigation fraction
            write(LIS_logunit,*) &
             "[INFO] Noah-MP.5.0 reading parameter IRFRACT from ", &
             trim(LIS_rc%paramfile(n))
            call LIS_read_param(n, trim(NoahMP50_struc(n)%LDT_ncvar_irfract), placeholder)
            do t = 1, LIS_rc%npatch(n, mtype)
                col = LIS_surface(n, mtype)%tile(t)%col
                row = LIS_surface(n, mtype)%tile(t)%row
                NoahMP50_struc(n)%noahmp50(t)%irfract = placeholder(col, row)
            enddo

            ! read: sprinkler irrigation fraction
            write(LIS_logunit,*) &
             "[INFO] Noah-MP.5.0 reading parameter SIFRACT from ", &
             trim(LIS_rc%paramfile(n))
            call LIS_read_param(n, trim(NoahMP50_struc(n)%LDT_ncvar_sifract), placeholder)
            do t = 1, LIS_rc%npatch(n, mtype)
                col = LIS_surface(n, mtype)%tile(t)%col
                row = LIS_surface(n, mtype)%tile(t)%row
                NoahMP50_struc(n)%noahmp50(t)%sifract = placeholder(col, row)
            enddo

            ! read: micro/drip irrigation fraction
            write(LIS_logunit,*) &
             "[INFO] Noah-MP.5.0 reading parameter MIFRACT from ", &
             trim(LIS_rc%paramfile(n))
            call LIS_read_param(n, trim(NoahMP50_struc(n)%LDT_ncvar_mifract), placeholder)
            do t = 1, LIS_rc%npatch(n, mtype)
                col = LIS_surface(n, mtype)%tile(t)%col
                row = LIS_surface(n, mtype)%tile(t)%row
                NoahMP50_struc(n)%noahmp50(t)%mifract = placeholder(col, row)
            enddo

            ! read: flood irrigation fraction
            write(LIS_logunit,*) &
             "[INFO] Noah-MP.5.0 reading parameter FIFRACT from ", &
             trim(LIS_rc%paramfile(n))
            call LIS_read_param(n, trim(NoahMP50_struc(n)%LDT_ncvar_fifract), placeholder)
            do t = 1, LIS_rc%npatch(n, mtype)
                col = LIS_surface(n, mtype)%tile(t)%col
                row = LIS_surface(n, mtype)%tile(t)%row
                NoahMP50_struc(n)%noahmp50(t)%fifract = placeholder(col, row)
            enddo
        endif

        ! for tile drainage
        if(NoahMP50_struc(n)%tdrn_opt > 0) then
            ! read: tile drainage fraction
            write(LIS_logunit,*) &
             "[INFO] Noah-MP.5.0 reading parameter TDFRACT from ", &
             trim(LIS_rc%paramfile(n))
            call LIS_read_param(n, trim(NoahMP50_struc(n)%LDT_ncvar_tdfract), placeholder)
            do t = 1, LIS_rc%npatch(n, mtype)
                col = LIS_surface(n, mtype)%tile(t)%col
                row = LIS_surface(n, mtype)%tile(t)%row
                NoahMP50_struc(n)%noahmp50(t)%tdfract = placeholder(col, row)
            enddo
        endif
 
        if(NoahMP50_struc(n)%soil_opt .eq. 2) then 
            ! read: soilcL1
            write(LIS_logunit,*) "[INFO] Noah-MP.5.0 reading parameter SOILCL1 from ", &
                 trim(LIS_rc%paramfile(n))
            call LIS_read_param(n, trim(NoahMP50_struc(n)%LDT_ncvar_soilcL1), placeholder)
            do t = 1, LIS_rc%npatch(n, mtype)
                col = LIS_surface(n, mtype)%tile(t)%col
                row = LIS_surface(n, mtype)%tile(t)%row
                NoahMP50_struc(n)%noahmp50(t)%soilcl1 = placeholder(col, row)
            enddo 

            ! read: soilcL2
            write(LIS_logunit,*) "[INFO] Noah-MP.5.0 reading parameter SOILCL2 from ", &
                 trim(LIS_rc%paramfile(n))
            call LIS_read_param(n, trim(NoahMP50_struc(n)%LDT_ncvar_soilcL2), placeholder)
            do t = 1, LIS_rc%npatch(n, mtype)
                col = LIS_surface(n, mtype)%tile(t)%col
                row = LIS_surface(n, mtype)%tile(t)%row
                NoahMP50_struc(n)%noahmp50(t)%soilcl2 = placeholder(col, row)
            enddo 

            ! read: soilcL3
            write(LIS_logunit,*) "[INFO] Noah-MP.5.0 reading parameter SOILCL3 from ", &
                 trim(LIS_rc%paramfile(n))
            call LIS_read_param(n, trim(NoahMP50_struc(n)%LDT_ncvar_soilcL3), placeholder)
            do t = 1, LIS_rc%npatch(n, mtype)
                col = LIS_surface(n, mtype)%tile(t)%col
                row = LIS_surface(n, mtype)%tile(t)%row
                NoahMP50_struc(n)%noahmp50(t)%soilcl3 = placeholder(col, row)
            enddo 

            ! read: soilcL4
            write(LIS_logunit,*) "[INFO] Noah-MP.5.0 reading parameter SOILCL4 from ", &
                 trim(LIS_rc%paramfile(n))
            call LIS_read_param(n, trim(NoahMP50_struc(n)%LDT_ncvar_soilcL4), placeholder)
            do t = 1, LIS_rc%npatch(n, mtype)
                col = LIS_surface(n, mtype)%tile(t)%col
                row = LIS_surface(n, mtype)%tile(t)%row
                NoahMP50_struc(n)%noahmp50(t)%soilcl4 = placeholder(col, row)
            enddo 
        endif 

        ! for MMF groundwater
        if(NoahMP50_struc(n)%runsub_opt == 5) then
            ! read: efolding depth for transmissivity (m)
            write(LIS_logunit,*) &
             "[INFO] Noah-MP.5.0 reading parameter FDEPTH from ", &
             trim(LIS_rc%paramfile(n))
            call LIS_read_param(n, trim(NoahMP50_struc(n)%LDT_ncvar_fdepth), placeholder)
            do t = 1, LIS_rc%npatch(n, mtype)
                col = LIS_surface(n, mtype)%tile(t)%col
                row = LIS_surface(n, mtype)%tile(t)%row
                NoahMP50_struc(n)%noahmp50(t)%fdepth = placeholder(col, row)
            enddo

            ! read: equilibrium water table depth (m)
            write(LIS_logunit,*) &
             "[INFO] Noah-MP.5.0 reading parameter EQZWT from ", &
             trim(LIS_rc%paramfile(n))
            call LIS_read_param(n, trim(NoahMP50_struc(n)%LDT_ncvar_eqzwt), placeholder)
            do t = 1, LIS_rc%npatch(n, mtype)
                col = LIS_surface(n, mtype)%tile(t)%col
                row = LIS_surface(n, mtype)%tile(t)%row
                NoahMP50_struc(n)%noahmp50(t)%eqzwt = placeholder(col, row)
            enddo

            ! read: riverbed depth (m)
            write(LIS_logunit,*) &
             "[INFO] Noah-MP.5.0 reading parameter RIVERBED from ", &
             trim(LIS_rc%paramfile(n))
            call LIS_read_param(n, trim(NoahMP50_struc(n)%LDT_ncvar_riverbed), placeholder)
            do t = 1, LIS_rc%npatch(n, mtype)
                col = LIS_surface(n, mtype)%tile(t)%col
                row = LIS_surface(n, mtype)%tile(t)%row
                NoahMP50_struc(n)%noahmp50(t)%riverbed = placeholder(col, row)
            enddo

            ! read: climatology recharge
            write(LIS_logunit,*) &
             "[INFO] Noah-MP.5.0 reading parameter RECHCLIM from ", &
             trim(LIS_rc%paramfile(n))
            call LIS_read_param(n, trim(NoahMP50_struc(n)%LDT_ncvar_rechclim), placeholder)
            do t = 1, LIS_rc%npatch(n, mtype)
                col = LIS_surface(n, mtype)%tile(t)%col
                row = LIS_surface(n, mtype)%tile(t)%row
                NoahMP50_struc(n)%noahmp50(t)%rechclim = placeholder(col, row)
            enddo
        endif

        !----------------------------------------------!
        ! MULTILEVEL reading spatial spatial parameters !
        !----------------------------------------------!
        write(LIS_logunit,*) "[INFO] Noah-MP.5.0 reading parameter SHDFAC_MONTHLY from ",&
             trim(LIS_rc%paramfile(n))
        do k = 1, 12
            call NoahMP50_read_MULTILEVEL_param(n, NoahMP50_struc(n)%LDT_ncvar_shdfac_monthly, k, placeholder)
            do t = 1, LIS_rc%npatch(n, mtype)
                col = LIS_surface(n, mtype)%tile(t)%col
                row = LIS_surface(n, mtype)%tile(t)%row
                NoahMP50_struc(n)%noahmp50(t)%shdfac_monthly(k) = placeholder(col, row)
            enddo 
        enddo 
        
        if(NoahMP50_struc(n)%soil_opt .eq. 3) then
            ! read: soilcomp
            write(LIS_logunit,*) "[INFO] Noah-MP.5.0 reading parameter SOILCOMP from ", &
                 trim(LIS_rc%paramfile(n))
            do k = 1, 8
                call NoahMP50_read_MULTILEVEL_param(n, NoahMP50_struc(n)%LDT_ncvar_soilcomp, k, placeholder)
                do t = 1, LIS_rc%npatch(n, mtype)
                    col = LIS_surface(n, mtype)%tile(t)%col
                    row = LIS_surface(n, mtype)%tile(t)%row
                    NoahMP50_struc(n)%noahmp50(t)%soilcomp(k) = placeholder(col, row)
                enddo 
            enddo 
        endif
        deallocate(placeholder)

        !!!! read Noah-MP parameter tables 
        write(LIS_logunit,*) "[INFO] Noah-MP.5.0 parameter table (veg, soil, general): ", &
             trim(NoahMP50_struc(n)%noahmp_tbl_name)     
        write(LIS_logunit,*) "[INFO] Noah-MP.5.0 Landuse classification scheme: ", &
             trim(NoahMP50_struc(n)%landuse_scheme_name) 
        write(LIS_logunit,*) "[INFO] Noah-MP.5.0 Soil classification scheme: ",  &
             "STAS (default, cannot change)" 
        call NoahmpReadTable(trim(NoahMP50_struc(n)%landuse_scheme_name), &
                             trim(NoahMP50_struc(n)%noahmp_tbl_name))

        do t=1,LIS_rc%npatch(n,mtype)
           SOILTYP = NoahMP50_struc(n)%noahmp50(t)%soiltype
           VEGTYP  = NoahMP50_struc(n)%noahmp50(t)%vegetype
           SLOPETYP     = 1          ! set underground runoff slope term
           SOILCOLOR    = 4          ! soil color: assuming a middle color category ?????????      
           CROPTYPE     = 0 
           if (NoahMP50_struc(n)%crop_opt > 0 .and. VEGTYP == NoahmpIO%ISCROP_TABLE) &
               CROPTYPE = NoahmpIO%DEFAULT_CROP_TABLE
           call TRANSFER_MP_PARAMETERS_NEW(VEGTYP,SOILTYP,SLOPETYP,SOILCOLOR,CROPTYPE,&
                NoahMP50_struc(n)%noahmp50(t)%param)
        enddo
   
        ! optional read of Optimized parameters
        call NoahMP50_read_OPT_parameters()

        !-------- initialize NoahmpIO 1-D interface variables
        NoahmpIO%xstart = 1
        NoahmpIO%xend   = 1
        NoahmpIO%ystart = 1
        NoahmpIO%yend   = 1
        NoahmpIO%ids    = NoahmpIO%xstart
        NoahmpIO%ide    = NoahmpIO%xend
        NoahmpIO%jds    = NoahmpIO%ystart
        NoahmpIO%jde    = NoahmpIO%yend
        NoahmpIO%kds    = 1
        NoahmpIO%kde    = 2
        NoahmpIO%its    = NoahmpIO%xstart
        NoahmpIO%ite    = NoahmpIO%xend
        NoahmpIO%jts    = NoahmpIO%ystart
        NoahmpIO%jte    = NoahmpIO%yend
        NoahmpIO%kts    = 1
        NoahmpIO%kte    = 2
        NoahmpIO%ims    = NoahmpIO%xstart
        NoahmpIO%ime    = NoahmpIO%xend
        NoahmpIO%jms    = NoahmpIO%ystart
        NoahmpIO%jme    = NoahmpIO%yend
        NoahmpIO%kms    = 1
        NoahmpIO%kme    = 2
        NoahmpIO%nsoil  = NoahMP50_struc(n)%nsoil
        NoahmpIO%nsnow  = NoahMP50_struc(n)%nsnow

        call NoahmpIOVarInitDefault(NoahmpIO) ! initialize NoahmpIO to undefined/default value
        !-------- NoahmpIO init complete

     enddo

end subroutine NoahMP50_setup

!BOP
!
! !ROUTINE: NoahMP50_read_MULTILEVEL_param
!  \label{read_MULTILEVEL_param}
!
! !REVISION HISTORY:
!  05/01/23: Cenlin He; update to work with refactored Noah-MP (v5.0 and later)
!
! !INTERFACE:
subroutine NoahMP50_read_MULTILEVEL_param(n, ncvar_name, level, placeholder)
! !USES:
#if (defined USE_NETCDF3 || defined USE_NETCDF4)
    use netcdf
#endif
    use LIS_coreMod, only : LIS_rc, LIS_localPet,   &   
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


    integer       :: ios, nid, param_ID, nc_ID, nr_ID, dimids(3)
    integer       :: nc, nr, nlevel
    real, pointer :: level_data(:, :, :)
    logical       :: file_exists

#if (defined USE_NETCDF3)
  write(LIS_logunit,*) "[ERR] NoahMP50_read_MULTILEVEL_param requires NetCDF4"
  call LIS_endrun()
#endif

#if (defined USE_NETCDF4)
    inquire(file=LIS_rc%paramfile(n), exist=file_exists)
    if(file_exists) then
        write(LIS_logunit, *) '[INFO] Reading '//trim(ncvar_name)//' map for level ', level

        ! open NetCDF parameter file
        ios = nf90_open(path=trim(LIS_rc%paramfile(n)), mode=NF90_NOWRITE, ncid=nid)
        call LIS_verify(ios, 'Error in nf90_open in NoahMP50_read_MULTILEVEL_param')

        ! inquire the ID of east-west dimension
        ios = nf90_inq_dimid(nid, 'east_west', nc_ID)
        call LIS_verify(ios, 'Error in nf90_inq_dimid in NoahMP50_read_MULTILEVEL_param')

        ! inquire the ID of north-south dimension
        ios = nf90_inq_dimid(nid, 'north_south', nr_ID)
        call LIS_verify(ios, 'Error in nf90_inq_dimid in NoahMP50_read_MULTILEVEL_param')

        ! inquire the length of east-west dimension
        ios = nf90_inquire_dimension(nid, nc_ID, len=nc)
        call LIS_verify(ios, 'Error in nf90_inquire_dimension in NoahMP50_read_MULTILEVEL_param')

        ! inquire the length of north-south dimension
        ios = nf90_inquire_dimension(nid, nr_ID, len=nr)
        call LIS_verify(ios, 'Error in nf90_inquire_dimension in NoahMP50_read_MULTILEVEL_param')

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
        call LIS_verify(ios, 'Error in nf90_get_var in NoahMP50_read_MULTILEVEL_param')

        ! close netcdf file 
        ios = nf90_close(nid)
        call LIS_verify(ios, 'Error in nf90_close in NoahMP50_read_MULTILEVEL_param')

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
#endif

 end subroutine NoahMP50_read_MULTILEVEL_param
                                          
SUBROUTINE TRANSFER_MP_PARAMETERS_NEW(VEGTYPE,SOILTYPE,SLOPETYPE,SOILCOLOR,CROPTYPE,parameters)

  use NoahmpIOVarType
  use LisNoahmpParamType

  implicit none

  INTEGER, INTENT(IN)    :: VEGTYPE
  INTEGER, INTENT(IN)    :: SOILTYPE(4)
  INTEGER, INTENT(IN)    :: SLOPETYPE
  INTEGER, INTENT(IN)    :: SOILCOLOR
  INTEGER, INTENT(IN)    :: CROPTYPE
    
  type(LisNoahmpParam_type), intent(inout) :: parameters
    
  REAL    :: FRZFACT
  INTEGER :: ISOIL

  parameters%ISWATER   = NoahmpIO%ISWATER_TABLE
  parameters%ISBARREN  = NoahmpIO%ISBARREN_TABLE
  parameters%ISICE     = NoahmpIO%ISICE_TABLE
  parameters%ISCROP    = NoahmpIO%ISCROP_TABLE
  parameters%EBLFOREST = NoahmpIO%EBLFOREST_TABLE

  parameters%URBAN_FLAG = .FALSE.
  IF( VEGTYPE == NoahmpIO%ISURBAN_TABLE .or. VEGTYPE >= NoahmpIO%URBTYPE_beg ) THEN
     parameters%URBAN_FLAG = .TRUE.
  ENDIF

!------------------------------------------------------------------------------------------!
! Transfer veg parameters
!------------------------------------------------------------------------------------------!

  parameters%CH2OP   =  NoahmpIO%CH2OP_TABLE(VEGTYPE)       !maximum intercepted h2o per unit lai+sai (mm)
  parameters%DLEAF   =  NoahmpIO%DLEAF_TABLE(VEGTYPE)       !characteristic leaf dimension (m)
  parameters%Z0MVT   =  NoahmpIO%Z0MVT_TABLE(VEGTYPE)       !momentum roughness length (m)
  parameters%HVT     =    NoahmpIO%HVT_TABLE(VEGTYPE)       !top of canopy (m)
  parameters%HVB     =    NoahmpIO%HVB_TABLE(VEGTYPE)       !bottom of canopy (m)
  parameters%DEN     =    NoahmpIO%DEN_TABLE(VEGTYPE)       !tree density (no. of trunks per m2)
  parameters%RC      =     NoahmpIO%RC_TABLE(VEGTYPE)       !tree crown radius (m)
  parameters%MFSNO   =  NoahmpIO%MFSNO_TABLE(VEGTYPE)       !snowmelt m parameter ()
  parameters%SCFFAC  = NoahmpIO%SCFFAC_TABLE(VEGTYPE)       !snow cover factor (m) (replace original hard-coded 2.5*z0 in SCF formulation)
  parameters%CBIOM   =  NoahmpIO%CBIOM_TABLE(VEGTYPE)       !canopy biomass heat capacity parameter (m)
  parameters%SAIM    =   NoahmpIO%SAIM_TABLE(VEGTYPE,:)     !monthly stem area index, one-sided
  parameters%LAIM    =   NoahmpIO%LAIM_TABLE(VEGTYPE,:)     !monthly leaf area index, one-sided
  parameters%SLA     =    NoahmpIO%SLA_TABLE(VEGTYPE)       !single-side leaf area per Kg [m2/kg]
  parameters%DILEFC  = NoahmpIO%DILEFC_TABLE(VEGTYPE)       !coeficient for leaf stress death [1/s]
  parameters%DILEFW  = NoahmpIO%DILEFW_TABLE(VEGTYPE)       !coeficient for leaf stress death [1/s]
  parameters%FRAGR   =  NoahmpIO%FRAGR_TABLE(VEGTYPE)       !fraction of growth respiration  !original was 0.3 
  parameters%LTOVRC  = NoahmpIO%LTOVRC_TABLE(VEGTYPE)       !leaf turnover [1/s]
  parameters%C3PSN   =  NoahmpIO%C3PSN_TABLE(VEGTYPE)       !photosynthetic pathway: 0. = c4, 1. = c3
  parameters%KC25    =   NoahmpIO%KC25_TABLE(VEGTYPE)       !co2 michaelis-menten constant at 25c (pa)
  parameters%AKC     =    NoahmpIO%AKC_TABLE(VEGTYPE)       !q10 for kc25
  parameters%KO25    =   NoahmpIO%KO25_TABLE(VEGTYPE)       !o2 michaelis-menten constant at 25c (pa)
  parameters%AKO     =    NoahmpIO%AKO_TABLE(VEGTYPE)       !q10 for ko25
  parameters%VCMX25  = NoahmpIO%VCMX25_TABLE(VEGTYPE)       !maximum rate of carboxylation at 25c (umol co2/m**2/s)
  parameters%AVCMX   =  NoahmpIO%AVCMX_TABLE(VEGTYPE)       !q10 for vcmx25
  parameters%BP      =     NoahmpIO%BP_TABLE(VEGTYPE)       !minimum leaf conductance (umol/m**2/s)
  parameters%MP      =     NoahmpIO%MP_TABLE(VEGTYPE)       !slope of conductance-to-photosynthesis relationship
  parameters%QE25    =   NoahmpIO%QE25_TABLE(VEGTYPE)       !quantum efficiency at 25c (umol co2 / umol photon)
  parameters%AQE     =    NoahmpIO%AQE_TABLE(VEGTYPE)       !q10 for qe25
  parameters%RMF25   =  NoahmpIO%RMF25_TABLE(VEGTYPE)       !leaf maintenance respiration at 25c (umol co2/m**2/s)
  parameters%RMS25   =  NoahmpIO%RMS25_TABLE(VEGTYPE)       !stem maintenance respiration at 25c (umol co2/kg bio/s)
  parameters%RMR25   =  NoahmpIO%RMR25_TABLE(VEGTYPE)       !root maintenance respiration at 25c (umol co2/kg bio/s)
  parameters%ARM     =    NoahmpIO%ARM_TABLE(VEGTYPE)       !q10 for maintenance respiration
  parameters%FOLNMX  = NoahmpIO%FOLNMX_TABLE(VEGTYPE)       !foliage nitrogen concentration when f(n)=1 (%)
  parameters%TMIN    =   NoahmpIO%TMIN_TABLE(VEGTYPE)       !minimum temperature for photosynthesis (k)
  parameters%XL      =     NoahmpIO%XL_TABLE(VEGTYPE)       !leaf/stem orientation index
  parameters%RHOL    =   NoahmpIO%RHOL_TABLE(VEGTYPE,:)     !leaf reflectance: 1=vis, 2=nir
  parameters%RHOS    =   NoahmpIO%RHOS_TABLE(VEGTYPE,:)     !stem reflectance: 1=vis, 2=nir
  parameters%TAUL    =   NoahmpIO%TAUL_TABLE(VEGTYPE,:)     !leaf transmittance: 1=vis, 2=nir
  parameters%TAUS    =   NoahmpIO%TAUS_TABLE(VEGTYPE,:)     !stem transmittance: 1=vis, 2=nir
  parameters%MRP     =    NoahmpIO%MRP_TABLE(VEGTYPE)       !microbial respiration parameter (umol co2 /kg c/ s)
  parameters%CWPVT   =  NoahmpIO%CWPVT_TABLE(VEGTYPE)       !empirical canopy wind parameter
  parameters%WRRAT   =  NoahmpIO%WRRAT_TABLE(VEGTYPE)       !wood to non-wood ratio
  parameters%WDPOOL  = NoahmpIO%WDPOOL_TABLE(VEGTYPE)       !wood pool (switch 1 or 0) depending on woody or not [-]
  parameters%TDLEF   =  NoahmpIO%TDLEF_TABLE(VEGTYPE)       !characteristic T for leaf freezing [K]
  parameters%NROOT   =  NoahmpIO%NROOT_TABLE(VEGTYPE)       !number of soil layers with root present
  parameters%RGL     =    NoahmpIO%RGL_TABLE(VEGTYPE)       !Parameter used in radiation stress function
  parameters%RSMIN   =     NoahmpIO%RS_TABLE(VEGTYPE)       !Minimum stomatal resistance [s m-1]
  parameters%HS      =     NoahmpIO%HS_TABLE(VEGTYPE)       !Parameter used in vapor pressure deficit function
  parameters%TOPT    =   NoahmpIO%TOPT_TABLE(VEGTYPE)       !Optimum transpiration air temperature [K]
  parameters%RSMAX   =  NoahmpIO%RSMAX_TABLE(VEGTYPE)       !Maximal stomatal resistance [s m-1]
  parameters%RTOVRC  = NoahmpIO%RTOVRC_TABLE(VEGTYPE)       !root turnover coefficient [1/s]
  parameters%RSWOODC = NoahmpIO%RSWOODC_TABLE(VEGTYPE)     !wood respiration coeficient [1/s]
  parameters%BF      =     NoahmpIO%BF_TABLE(VEGTYPE)       !parameter for present wood allocation [-]
  parameters%WSTRC   =  NoahmpIO%WSTRC_TABLE(VEGTYPE)       !water stress coeficient [-]
  parameters%LAIMIN  = NoahmpIO%LAIMIN_TABLE(VEGTYPE)       !minimum leaf area index [m2/m2]
  parameters%XSAMIN  = NoahmpIO%XSAMIN_TABLE(VEGTYPE)       !minimum stem area index [m2/m2]

!------------------------------------------------------------------------------------------!
! Transfer rad parameters
!------------------------------------------------------------------------------------------!

   parameters%ALBSAT = NoahmpIO%ALBSAT_TABLE(SOILCOLOR,:)
   parameters%ALBDRY = NoahmpIO%ALBDRY_TABLE(SOILCOLOR,:)
   parameters%ALBICE = NoahmpIO%ALBICE_TABLE
   parameters%ALBLAK = NoahmpIO%ALBLAK_TABLE               
   parameters%OMEGAS = NoahmpIO%OMEGAS_TABLE
   parameters%BETADS = NoahmpIO%BETADS_TABLE
   parameters%BETAIS = NoahmpIO%BETAIS_TABLE
   parameters%EG     = NoahmpIO%EG_TABLE
   parameters%EICE   = NoahmpIO%EICE_TABLE

!------------------------------------------------------------------------------------------!
! Transfer crop parameters
!------------------------------------------------------------------------------------------!

  IF(CROPTYPE > 0) THEN
   parameters%PLTDAY    =    NoahmpIO%PLTDAY_TABLE(CROPTYPE)    ! Planting date
   parameters%HSDAY     =     NoahmpIO%HSDAY_TABLE(CROPTYPE)    ! Harvest date
   parameters%PLANTPOP  =  NoahmpIO%PLANTPOP_TABLE(CROPTYPE)    ! Plant density [per ha] - used?
   parameters%IRRI      =      NoahmpIO%IRRI_TABLE(CROPTYPE)    ! Irrigation strategy 0= non-irrigation 1=irrigation (no water-stress)
   parameters%GDDTBASE  =  NoahmpIO%GDDTBASE_TABLE(CROPTYPE)    ! Base temperature for GDD accumulation [C]
   parameters%GDDTCUT   =   NoahmpIO%GDDTCUT_TABLE(CROPTYPE)    ! Upper temperature for GDD accumulation [C]
   parameters%GDDS1     =     NoahmpIO%GDDS1_TABLE(CROPTYPE)    ! GDD from seeding to emergence
   parameters%GDDS2     =     NoahmpIO%GDDS2_TABLE(CROPTYPE)    ! GDD from seeding to initial vegetative 
   parameters%GDDS3     =     NoahmpIO%GDDS3_TABLE(CROPTYPE)    ! GDD from seeding to post vegetative 
   parameters%GDDS4     =     NoahmpIO%GDDS4_TABLE(CROPTYPE)    ! GDD from seeding to intial reproductive
   parameters%GDDS5     =     NoahmpIO%GDDS5_TABLE(CROPTYPE)    ! GDD from seeding to pysical maturity 
   parameters%C3PSN     =    NoahmpIO%C3PSNI_TABLE(CROPTYPE)
   parameters%KC25      =     NoahmpIO%KC25I_TABLE(CROPTYPE)
   parameters%AKC       =      NoahmpIO%AKCI_TABLE(CROPTYPE)
   parameters%KO25      =     NoahmpIO%KO25I_TABLE(CROPTYPE)
   parameters%AKO       =      NoahmpIO%AKOI_TABLE(CROPTYPE)
   parameters%AVCMX     =    NoahmpIO%AVCMXI_TABLE(CROPTYPE)
   parameters%VCMX25    =   NoahmpIO%VCMX25I_TABLE(CROPTYPE)
   parameters%BP        =       NoahmpIO%BPI_TABLE(CROPTYPE)
   parameters%MP        =       NoahmpIO%MPI_TABLE(CROPTYPE)
   parameters%FOLNMX    =   NoahmpIO%FOLNMXI_TABLE(CROPTYPE)
   parameters%QE25      =     NoahmpIO%QE25I_TABLE(CROPTYPE)   
   parameters%AREF      =      NoahmpIO%AREF_TABLE(CROPTYPE)    ! reference maximum CO2 assimulation rate 
   parameters%PSNRF     =     NoahmpIO%PSNRF_TABLE(CROPTYPE)    ! CO2 assimulation reduction factor(0-1) (caused by non-modeling part,e.g.pest,weeds)
   parameters%I2PAR     =     NoahmpIO%I2PAR_TABLE(CROPTYPE)    ! Fraction of incoming solar radiation to photosynthetically active radiation
   parameters%TASSIM0   =   NoahmpIO%TASSIM0_TABLE(CROPTYPE)    ! Minimum temperature for CO2 assimulation [C]
   parameters%TASSIM1   =   NoahmpIO%TASSIM1_TABLE(CROPTYPE)    ! CO2 assimulation linearly increasing until temperature reaches T1 [C]
   parameters%TASSIM2   =   NoahmpIO%TASSIM2_TABLE(CROPTYPE)    ! CO2 assmilation rate remain at Aref until temperature reaches T2 [C]
   parameters%K         =         NoahmpIO%K_TABLE(CROPTYPE)    ! light extinction coefficient
   parameters%EPSI      =      NoahmpIO%EPSI_TABLE(CROPTYPE)    ! initial light use efficiency
   parameters%Q10MR     =     NoahmpIO%Q10MR_TABLE(CROPTYPE)    ! q10 for maintainance respiration
   parameters%LEFREEZ   =   NoahmpIO%LEFREEZ_TABLE(CROPTYPE)    ! characteristic T for leaf freezing [K]
   parameters%DILE_FC   =   NoahmpIO%DILE_FC_TABLE(CROPTYPE,:)  ! coeficient for temperature leaf stress death [1/s]
   parameters%DILE_FW   =   NoahmpIO%DILE_FW_TABLE(CROPTYPE,:)  ! coeficient for water leaf stress death [1/s]
   parameters%FRA_GR    =    NoahmpIO%FRA_GR_TABLE(CROPTYPE)    ! fraction of growth respiration
   parameters%LF_OVRC   =   NoahmpIO%LF_OVRC_TABLE(CROPTYPE,:)  ! fraction of leaf turnover  [1/s]
   parameters%ST_OVRC   =   NoahmpIO%ST_OVRC_TABLE(CROPTYPE,:)  ! fraction of stem turnover  [1/s]
   parameters%RT_OVRC   =   NoahmpIO%RT_OVRC_TABLE(CROPTYPE,:)  ! fraction of root tunrover  [1/s]
   parameters%LFMR25    =    NoahmpIO%LFMR25_TABLE(CROPTYPE)    ! leaf maintenance respiration at 25C [umol CO2/m**2  /s]
   parameters%STMR25    =    NoahmpIO%STMR25_TABLE(CROPTYPE)    ! stem maintenance respiration at 25C [umol CO2/kg bio/s]
   parameters%RTMR25    =    NoahmpIO%RTMR25_TABLE(CROPTYPE)    ! root maintenance respiration at 25C [umol CO2/kg bio/s]
   parameters%GRAINMR25 = NoahmpIO%GRAINMR25_TABLE(CROPTYPE)    ! grain maintenance respiration at 25C [umol CO2/kg bio/s]
   parameters%LFPT      =      NoahmpIO%LFPT_TABLE(CROPTYPE,:)  ! fraction of carbohydrate flux to leaf
   parameters%STPT      =      NoahmpIO%STPT_TABLE(CROPTYPE,:)  ! fraction of carbohydrate flux to stem
   parameters%RTPT      =      NoahmpIO%RTPT_TABLE(CROPTYPE,:)  ! fraction of carbohydrate flux to root
   parameters%GRAINPT   =   NoahmpIO%GRAINPT_TABLE(CROPTYPE,:)  ! fraction of carbohydrate flux to grain
   parameters%LFCT      =      NoahmpIO%LFCT_TABLE(CROPTYPE,:)  ! fraction of carbohydrate translocation from leaf to grain
   parameters%STCT      =      NoahmpIO%STCT_TABLE(CROPTYPE,:)  ! fraction of carbohydrate translocation from stem to grain
   parameters%RTCT      =      NoahmpIO%RTCT_TABLE(CROPTYPE,:)  ! fraction of carbohydrate translocation from root to grain
   parameters%BIO2LAI   =   NoahmpIO%BIO2LAI_TABLE(CROPTYPE)    ! leaf are per living leaf biomass [m^2/kg]
  END IF

!------------------------------------------------------------------------------------------!
! Transfer global parameters
!------------------------------------------------------------------------------------------!

   parameters%CO2              =              NoahmpIO%CO2_TABLE
   parameters%O2               =               NoahmpIO%O2_TABLE
   parameters%TIMEAN           =           NoahmpIO%TIMEAN_TABLE
   parameters%FSATMX           =           NoahmpIO%FSATMX_TABLE
   parameters%Z0SNO            =            NoahmpIO%Z0SNO_TABLE
   parameters%SSI              =              NoahmpIO%SSI_TABLE
   parameters%SNOW_RET_FAC     =     NoahmpIO%SNOW_RET_FAC_TABLE
   parameters%SNOW_EMIS        =        NoahmpIO%SNOW_EMIS_TABLE
   parameters%SWEMX            =            NoahmpIO%SWEMX_TABLE
   parameters%RSURF_SNOW       =       NoahmpIO%RSURF_SNOW_TABLE
   parameters%TAU0             =             NoahmpIO%TAU0_TABLE
   parameters%GRAIN_GROWTH     =     NoahmpIO%GRAIN_GROWTH_TABLE
   parameters%EXTRA_GROWTH     =     NoahmpIO%EXTRA_GROWTH_TABLE
   parameters%DIRT_SOOT        =        NoahmpIO%DIRT_SOOT_TABLE
   parameters%BATS_COSZ        =        NoahmpIO%BATS_COSZ_TABLE
   parameters%BATS_VIS_NEW     =     NoahmpIO%BATS_VIS_NEW_TABLE
   parameters%BATS_NIR_NEW     =     NoahmpIO%BATS_NIR_NEW_TABLE
   parameters%BATS_VIS_AGE     =     NoahmpIO%BATS_VIS_AGE_TABLE
   parameters%BATS_NIR_AGE     =     NoahmpIO%BATS_NIR_AGE_TABLE
   parameters%BATS_VIS_DIR     =     NoahmpIO%BATS_VIS_DIR_TABLE
   parameters%BATS_NIR_DIR     =     NoahmpIO%BATS_NIR_DIR_TABLE
   parameters%RSURF_EXP        =        NoahmpIO%RSURF_EXP_TABLE
   parameters%C2_SNOWCOMPACT   =   NoahmpIO%C2_SNOWCOMPACT_TABLE
   parameters%C3_SNOWCOMPACT   =   NoahmpIO%C3_SNOWCOMPACT_TABLE 
   parameters%C4_SNOWCOMPACT   =   NoahmpIO%C4_SNOWCOMPACT_TABLE 
   parameters%C5_SNOWCOMPACT   =   NoahmpIO%C5_SNOWCOMPACT_TABLE
   parameters%DM_SNOWCOMPACT   =   NoahmpIO%DM_SNOWCOMPACT_TABLE
   parameters%ETA0_SNOWCOMPACT = NoahmpIO%ETA0_SNOWCOMPACT_TABLE
   parameters%SNLIQMAXFRAC     =     NoahmpIO%SNLIQMAXFRAC_TABLE
   parameters%SWEMAXGLA        =        NoahmpIO%SWEMAXGLA_TABLE
   parameters%WSLMAX           =           NoahmpIO%WSLMAX_TABLE
   parameters%ROUS             =             NoahmpIO%ROUS_TABLE
   parameters%CMIC             =             NoahmpIO%CMIC_TABLE
   parameters%SNOWDEN_MAX      =      NoahmpIO%SNOWDEN_MAX_TABLE
   parameters%CLASS_ALB_REF    =    NoahmpIO%CLASS_ALB_REF_TABLE
   parameters%CLASS_SNO_AGE    =    NoahmpIO%CLASS_SNO_AGE_TABLE
   parameters%CLASS_ALB_NEW    =    NoahmpIO%CLASS_ALB_NEW_TABLE
   parameters%PSIWLT           =           NoahmpIO%PSIWLT_TABLE
   parameters%Z0SOIL           =           NoahmpIO%Z0SOIL_TABLE
   parameters%Z0LAKE           =           NoahmpIO%Z0LAKE_TABLE

! ----------------------------------------------------------------------
!  Transfer irrigation parameters
! ----------------------------------------------------------------------
   parameters%IRR_HAR     =    NoahmpIO%IRR_HAR_TABLE
   parameters%IRR_FRAC    =   NoahmpIO%IRR_FRAC_TABLE
   parameters%IRR_LAI     =    NoahmpIO%IRR_LAI_TABLE
   parameters%IRR_MAD     =    NoahmpIO%IRR_MAD_TABLE 
   parameters%FILOSS      =     NoahmpIO%FILOSS_TABLE 
   parameters%SPRIR_RATE  = NoahmpIO%SPRIR_RATE_TABLE
   parameters%MICIR_RATE  = NoahmpIO%MICIR_RATE_TABLE
   parameters%FIRTFAC     =    NoahmpIO%FIRTFAC_TABLE
   parameters%IR_RAIN     =    NoahmpIO%IR_RAIN_TABLE

! ----------------------------------------------------------------------
!  Transfer tile drainage parameters
! ----------------------------------------------------------------------
   parameters%DRAIN_LAYER_OPT = NoahmpIO%DRAIN_LAYER_OPT_TABLE
   parameters%TD_DEPTH        =        NoahmpIO%TD_DEPTH_TABLE(SOILTYPE(1))
   parameters%TDSMC_FAC       =       NoahmpIO%TDSMC_FAC_TABLE(SOILTYPE(1))
   parameters%TD_DC           =           NoahmpIO%TD_DC_TABLE(SOILTYPE(1))
   parameters%TD_DCOEF        =        NoahmpIO%TD_DCOEF_TABLE(SOILTYPE(1))
   parameters%TD_D            =            NoahmpIO%TD_D_TABLE(SOILTYPE(1))
   parameters%TD_ADEPTH       =       NoahmpIO%TD_ADEPTH_TABLE(SOILTYPE(1))
   parameters%TD_RADI         =         NoahmpIO%TD_RADI_TABLE(SOILTYPE(1))
   parameters%TD_SPAC         =         NoahmpIO%TD_SPAC_TABLE(SOILTYPE(1))
   parameters%TD_DDRAIN       =       NoahmpIO%TD_DDRAIN_TABLE(SOILTYPE(1))
   parameters%KLAT_FAC        =        NoahmpIO%KLAT_FAC_TABLE(SOILTYPE(1))

! ----------------------------------------------------------------------
!  Transfer soil parameters
! ----------------------------------------------------------------------

    do isoil = 1, size(soiltype)
      parameters%BEXP(isoil)   = NoahmpIO%BEXP_TABLE   (SOILTYPE(isoil))
      parameters%DKSAT(isoil)  = NoahmpIO%DKSAT_TABLE  (SOILTYPE(isoil))
      parameters%DWSAT(isoil)  = NoahmpIO%DWSAT_TABLE  (SOILTYPE(isoil))
      parameters%PSISAT(isoil) = NoahmpIO%PSISAT_TABLE (SOILTYPE(isoil))
      parameters%QUARTZ(isoil) = NoahmpIO%QUARTZ_TABLE (SOILTYPE(isoil))
      parameters%SMCDRY(isoil) = NoahmpIO%SMCDRY_TABLE (SOILTYPE(isoil))
      parameters%SMCMAX(isoil) = NoahmpIO%SMCMAX_TABLE (SOILTYPE(isoil))
      parameters%SMCREF(isoil) = NoahmpIO%SMCREF_TABLE (SOILTYPE(isoil))
      parameters%SMCWLT(isoil) = NoahmpIO%SMCWLT_TABLE (SOILTYPE(isoil))
    end do
    
    parameters%BVIC   = NoahmpIO%BVIC_TABLE(SOILTYPE(1))
    parameters%AXAJ   = NoahmpIO%AXAJ_TABLE(SOILTYPE(1))
    parameters%BXAJ   = NoahmpIO%BXAJ_TABLE(SOILTYPE(1))
    parameters%XXAJ   = NoahmpIO%XXAJ_TABLE(SOILTYPE(1))
    parameters%BDVIC  = NoahmpIO%BDVIC_TABLE(SOILTYPE(1))
    parameters%GDVIC  = NoahmpIO%GDVIC_TABLE(SOILTYPE(1))
    parameters%BBVIC  = NoahmpIO%BBVIC_TABLE(SOILTYPE(1))

! ----------------------------------------------------------------------
! Transfer GENPARM parameters
! ----------------------------------------------------------------------
    parameters%CSOIL  = NoahmpIO%CSOIL_TABLE
    parameters%ZBOT   = NoahmpIO%ZBOT_TABLE
    parameters%CZIL   = NoahmpIO%CZIL_TABLE
    parameters%REFDK  = NoahmpIO%REFDK_TABLE
    parameters%REFKDT = NoahmpIO%REFKDT_TABLE
    parameters%FRZK   = NoahmpIO%FRZK_TABLE
    parameters%KDT    = parameters%REFKDT * parameters%DKSAT(1) / parameters%REFDK
    parameters%SLOPE  = NoahmpIO%SLOPE_TABLE(SLOPETYPE)

    IF(parameters%URBAN_FLAG)THEN  ! Hardcoding some urban parameters for soil
       parameters%SMCMAX = 0.45 
       parameters%SMCREF = 0.42 
       parameters%SMCWLT = 0.40 
       parameters%SMCDRY = 0.40 
       parameters%CSOIL  = 3.E6
    ENDIF

! adjust FRZK parameter to actual soil type: FRZK * FRZFACT
    IF(SOILTYPE(1) /= 14) then
      FRZFACT = (parameters%SMCMAX(1) / parameters%SMCREF(1)) * (0.412 / 0.468)
      parameters%FRZX = parameters%FRZK * FRZFACT
    END IF

    parameters%mxsnalb = 0.84
    parameters%mnsnalb = 0.55
    parameters%sndecayexp = 0.01
    parameters%t_ulimit = 2.5
    parameters%t_mlimit = 2.0
    parameters%t_llimit = 0.5
    parameters%snowf_scalef = 1.0
   
 END SUBROUTINE TRANSFER_MP_PARAMETERS_NEW
