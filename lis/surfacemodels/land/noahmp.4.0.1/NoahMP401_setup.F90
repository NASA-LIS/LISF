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
    integer       :: soilcolor, vegtyp, soiltyp(4), slopetyp, croptype    
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
             "[INFO] Noah-MP.4.0.1 reading parameter PLANTING from ", &
             trim(LIS_rc%paramfile(n))
            call LIS_read_param(n, trim(NOAHMP401_struc(n)%LDT_ncvar_planting), placeholder)
            do t = 1, LIS_rc%npatch(n, mtype)
                col = LIS_surface(n, mtype)%tile(t)%col
                row = LIS_surface(n, mtype)%tile(t)%row
                NOAHMP401_struc(n)%noahmp401(t)%planting = placeholder(col, row)
            enddo 

            ! read: harvest
            write(LIS_logunit,*) "[INFO] Noah-MP.4.0.1 reading parameter HARVEST from ",&
                 trim(LIS_rc%paramfile(n))
            call LIS_read_param(n, trim(NOAHMP401_struc(n)%LDT_ncvar_harvest), placeholder)
            do t = 1, LIS_rc%npatch(n, mtype)
                col = LIS_surface(n, mtype)%tile(t)%col
                row = LIS_surface(n, mtype)%tile(t)%row
                NOAHMP401_struc(n)%noahmp401(t)%harvest = placeholder(col, row)
            enddo 

            ! read: season_gdd
            write(LIS_logunit,*) "[INFO] Noah-MP.4.0.1 reading parameter SEASON_GDD from ", &
                 trim(LIS_rc%paramfile(n))
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
            write(LIS_logunit,*) "[INFO] Noah-MP.4.0.1 reading parameter SOILCL1 from ", &
                 trim(LIS_rc%paramfile(n))
            call LIS_read_param(n, trim(NOAHMP401_struc(n)%LDT_ncvar_soilcL1), placeholder)
            do t = 1, LIS_rc%npatch(n, mtype)
                col = LIS_surface(n, mtype)%tile(t)%col
                row = LIS_surface(n, mtype)%tile(t)%row
                NOAHMP401_struc(n)%noahmp401(t)%soilcl1 = placeholder(col, row)
            enddo 

            ! read: soilcL2
            write(LIS_logunit,*) "[INFO] Noah-MP.4.0.1 reading parameter SOILCL2 from ", &
                 trim(LIS_rc%paramfile(n))
            call LIS_read_param(n, trim(NOAHMP401_struc(n)%LDT_ncvar_soilcL2), placeholder)
            do t = 1, LIS_rc%npatch(n, mtype)
                col = LIS_surface(n, mtype)%tile(t)%col
                row = LIS_surface(n, mtype)%tile(t)%row
                NOAHMP401_struc(n)%noahmp401(t)%soilcl2 = placeholder(col, row)
            enddo 

            ! read: soilcL3
            write(LIS_logunit,*) "[INFO] Noah-MP.4.0.1 reading parameter SOILCL3 from ", &
                 trim(LIS_rc%paramfile(n))
            call LIS_read_param(n, trim(NOAHMP401_struc(n)%LDT_ncvar_soilcL3), placeholder)
            do t = 1, LIS_rc%npatch(n, mtype)
                col = LIS_surface(n, mtype)%tile(t)%col
                row = LIS_surface(n, mtype)%tile(t)%row
                NOAHMP401_struc(n)%noahmp401(t)%soilcl3 = placeholder(col, row)
            enddo 

            ! read: soilcL4
            write(LIS_logunit,*) "[INFO] Noah-MP.4.0.1 reading parameter SOILCL4 from ", &
                 trim(LIS_rc%paramfile(n))
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
        write(LIS_logunit,*) "[INFO] Noah-MP.4.0.1 reading parameter SHDFAC_MONTHLY from ",&
             trim(LIS_rc%paramfile(n))
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
            write(LIS_logunit,*) "[INFO] Noah-MP.4.0.1 reading parameter SOILCOMP from ", &
                 trim(LIS_rc%paramfile(n))
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
        write(LIS_logunit,*) "[INFO] Noah-MP.4.0.1 soil parameter table: ", &
             trim(NOAHMP401_struc(n)%soil_tbl_name)       
        write(LIS_logunit,*) "[INFO] Noah-MP.4.0.1 general parameter table: ",&
             trim(NOAHMP401_struc(n)%gen_tbl_name)        
        write(LIS_logunit,*) "[INFO] Noah-MP.4.0.1 MP parameter table: ",   &
             trim(NOAHMP401_struc(n)%noahmp_tbl_name)     
        write(LIS_logunit,*) "[INFO] Noah-MP.4.0.1 Landuse classification scheme: ", &
             trim(NOAHMP401_struc(n)%landuse_scheme_name) 
        write(LIS_logunit,*) "[INFO] Noah-MP.4.0.1 Soil classification scheme: ",  &
             "STAS (default, cannot change)" 
        call read_mp_veg_parameters(trim(NOAHMP401_struc(n)%landuse_scheme_name), &
             trim(NOAHMP401_struc(n)%noahmp_tbl_name))
        call read_mp_soil_parameters(trim(NOAHMP401_struc(n)%soil_tbl_name), &
             trim(NOAHMP401_struc(n)%gen_tbl_name))
        call read_mp_rad_parameters(trim(NOAHMP401_struc(n)%noahmp_tbl_name))
        call read_mp_global_parameters(trim(NOAHMP401_struc(n)%noahmp_tbl_name))
        call read_mp_crop_parameters(trim(NOAHMP401_struc(n)%noahmp_tbl_name))
        call read_mp_optional_parameters(trim(NOAHMP401_struc(n)%noahmp_tbl_name))

        do t=1,LIS_rc%npatch(n,mtype)
           soiltyp = NoahMP401_struc(n)%noahmp401(t)%soiltype
           vegtyp  = NoahMP401_struc(n)%noahmp401(t)%vegetype
           
           SLOPETYP     = 1          ! set underground runoff slope term
           SOILCOLOR    = 4          ! soil color: assuming a middle color category ?????????      
           CROPTYPE = 0 
           CALL TRANSFER_MP_PARAMETERS(VEGTYP,SOILTYP,SLOPETYP,SOILCOLOR,CROPTYPE,&
                NoahMP401_struc(n)%noahmp401(t)%param)

        enddo
   
     !optional read of Optimized parameters

        call NoahMP401_read_OPT_parameters()     
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
                                          
SUBROUTINE TRANSFER_MP_PARAMETERS(VEGTYPE,SOILTYPE,SLOPETYPE,SOILCOLOR,CROPTYPE,parameters)

  USE NOAHMP_TABLES_401
  USE MODULE_SF_NOAHMPLSM_401

  implicit none

  INTEGER, INTENT(IN)    :: VEGTYPE
  INTEGER, INTENT(IN)    :: SOILTYPE(4)
  INTEGER, INTENT(IN)    :: SLOPETYPE
  INTEGER, INTENT(IN)    :: SOILCOLOR
  INTEGER, INTENT(IN)    :: CROPTYPE
    
  type (noahmp_parameters), intent(inout) :: parameters
    
  REAL    :: REFDK
  REAL    :: REFKDT
  REAL    :: FRZK
  REAL    :: FRZFACT
  INTEGER :: ISOIL

  parameters%ISWATER   =   ISWATER_TABLE
  parameters%ISBARREN  =  ISBARREN_TABLE
  parameters%ISICE     =     ISICE_TABLE
  parameters%ISCROP    =    ISCROP_TABLE
  parameters%EBLFOREST = EBLFOREST_TABLE

  parameters%URBAN_FLAG = .FALSE.
  IF( VEGTYPE == ISURBAN_TABLE                  .or. VEGTYPE == LOW_DENSITY_RESIDENTIAL_TABLE  .or. &
      VEGTYPE == HIGH_DENSITY_RESIDENTIAL_TABLE .or. VEGTYPE == HIGH_INTENSITY_INDUSTRIAL_TABLE ) THEN
     parameters%URBAN_FLAG = .TRUE.
  ENDIF

!------------------------------------------------------------------------------------------!
! Transfer veg parameters
!------------------------------------------------------------------------------------------!

  parameters%CH2OP  =  CH2OP_TABLE(VEGTYPE)       !maximum intercepted h2o per unit lai+sai (mm)
  parameters%DLEAF  =  DLEAF_TABLE(VEGTYPE)       !characteristic leaf dimension (m)
  parameters%Z0MVT  =  Z0MVT_TABLE(VEGTYPE)       !momentum roughness length (m)
  parameters%HVT    =    HVT_TABLE(VEGTYPE)       !top of canopy (m)
  parameters%HVB    =    HVB_TABLE(VEGTYPE)       !bottom of canopy (m)
  parameters%DEN    =    DEN_TABLE(VEGTYPE)       !tree density (no. of trunks per m2)
  parameters%RC     =     RC_TABLE(VEGTYPE)       !tree crown radius (m)
  parameters%MFSNO  =  MFSNO_TABLE(VEGTYPE)       !snowmelt m parameter ()
  parameters%SAIM   =   SAIM_TABLE(VEGTYPE,:)     !monthly stem area index, one-sided
  parameters%LAIM   =   LAIM_TABLE(VEGTYPE,:)     !monthly leaf area index, one-sided
  parameters%SLA    =    SLA_TABLE(VEGTYPE)       !single-side leaf area per Kg [m2/kg]
  parameters%DILEFC = DILEFC_TABLE(VEGTYPE)       !coeficient for leaf stress death [1/s]
  parameters%DILEFW = DILEFW_TABLE(VEGTYPE)       !coeficient for leaf stress death [1/s]
  parameters%FRAGR  =  FRAGR_TABLE(VEGTYPE)       !fraction of growth respiration  !original was 0.3 
  parameters%LTOVRC = LTOVRC_TABLE(VEGTYPE)       !leaf turnover [1/s]

  parameters%C3PSN  =  C3PSN_TABLE(VEGTYPE)       !photosynthetic pathway: 0. = c4, 1. = c3
  parameters%KC25   =   KC25_TABLE(VEGTYPE)       !co2 michaelis-menten constant at 25c (pa)
  parameters%AKC    =    AKC_TABLE(VEGTYPE)       !q10 for kc25
  parameters%KO25   =   KO25_TABLE(VEGTYPE)       !o2 michaelis-menten constant at 25c (pa)
  parameters%AKO    =    AKO_TABLE(VEGTYPE)       !q10 for ko25
  parameters%VCMX25 = VCMX25_TABLE(VEGTYPE)       !maximum rate of carboxylation at 25c (umol co2/m**2/s)
  parameters%AVCMX  =  AVCMX_TABLE(VEGTYPE)       !q10 for vcmx25
  parameters%BP     =     BP_TABLE(VEGTYPE)       !minimum leaf conductance (umol/m**2/s)
  parameters%MP     =     MP_TABLE(VEGTYPE)       !slope of conductance-to-photosynthesis relationship
  parameters%QE25   =   QE25_TABLE(VEGTYPE)       !quantum efficiency at 25c (umol co2 / umol photon)
  parameters%AQE    =    AQE_TABLE(VEGTYPE)       !q10 for qe25
  parameters%RMF25  =  RMF25_TABLE(VEGTYPE)       !leaf maintenance respiration at 25c (umol co2/m**2/s)
  parameters%RMS25  =  RMS25_TABLE(VEGTYPE)       !stem maintenance respiration at 25c (umol co2/kg bio/s)
  parameters%RMR25  =  RMR25_TABLE(VEGTYPE)       !root maintenance respiration at 25c (umol co2/kg bio/s)
  parameters%ARM    =    ARM_TABLE(VEGTYPE)       !q10 for maintenance respiration
  parameters%FOLNMX = FOLNMX_TABLE(VEGTYPE)       !foliage nitrogen concentration when f(n)=1 (%)
  parameters%TMIN   =   TMIN_TABLE(VEGTYPE)       !minimum temperature for photosynthesis (k)

  parameters%XL     =     XL_TABLE(VEGTYPE)       !leaf/stem orientation index
  parameters%RHOL   =   RHOL_TABLE(VEGTYPE,:)     !leaf reflectance: 1=vis, 2=nir
  parameters%RHOS   =   RHOS_TABLE(VEGTYPE,:)     !stem reflectance: 1=vis, 2=nir
  parameters%TAUL   =   TAUL_TABLE(VEGTYPE,:)     !leaf transmittance: 1=vis, 2=nir
  parameters%TAUS   =   TAUS_TABLE(VEGTYPE,:)     !stem transmittance: 1=vis, 2=nir

  parameters%MRP    =    MRP_TABLE(VEGTYPE)       !microbial respiration parameter (umol co2 /kg c/ s)
  parameters%CWPVT  =  CWPVT_TABLE(VEGTYPE)       !empirical canopy wind parameter

  parameters%WRRAT  =  WRRAT_TABLE(VEGTYPE)       !wood to non-wood ratio
  parameters%WDPOOL = WDPOOL_TABLE(VEGTYPE)       !wood pool (switch 1 or 0) depending on woody or not [-]
  parameters%TDLEF  =  TDLEF_TABLE(VEGTYPE)       !characteristic T for leaf freezing [K]

  parameters%NROOT  =  NROOT_TABLE(VEGTYPE)       !number of soil layers with root present
  parameters%RGL    =    RGL_TABLE(VEGTYPE)       !Parameter used in radiation stress function
  parameters%RSMIN  =     RS_TABLE(VEGTYPE)       !Minimum stomatal resistance [s m-1]
  parameters%HS     =     HS_TABLE(VEGTYPE)       !Parameter used in vapor pressure deficit function
  parameters%TOPT   =   TOPT_TABLE(VEGTYPE)       !Optimum transpiration air temperature [K]
  parameters%RSMAX  =  RSMAX_TABLE(VEGTYPE)       !Maximal stomatal resistance [s m-1]

!------------------------------------------------------------------------------------------!
! Transfer rad parameters
!------------------------------------------------------------------------------------------!

   parameters%ALBSAT    = ALBSAT_TABLE(SOILCOLOR,:)
   parameters%ALBDRY    = ALBDRY_TABLE(SOILCOLOR,:)
   parameters%ALBICE    = ALBICE_TABLE
   parameters%ALBLAK    = ALBLAK_TABLE               
   parameters%OMEGAS    = OMEGAS_TABLE
   parameters%BETADS    = BETADS_TABLE
   parameters%BETAIS    = BETAIS_TABLE
   parameters%EG        = EG_TABLE

!------------------------------------------------------------------------------------------!
! Transfer crop parameters
!------------------------------------------------------------------------------------------!

  IF(CROPTYPE > 0) THEN
   parameters%PLTDAY    =    PLTDAY_TABLE(CROPTYPE)    ! Planting date
   parameters%HSDAY     =     HSDAY_TABLE(CROPTYPE)    ! Harvest date
   parameters%PLANTPOP  =  PLANTPOP_TABLE(CROPTYPE)    ! Plant density [per ha] - used?
   parameters%IRRI      =      IRRI_TABLE(CROPTYPE)    ! Irrigation strategy 0= non-irrigation 1=irrigation (no water-stress)
   parameters%GDDTBASE  =  GDDTBASE_TABLE(CROPTYPE)    ! Base temperature for GDD accumulation [C]
   parameters%GDDTCUT   =   GDDTCUT_TABLE(CROPTYPE)    ! Upper temperature for GDD accumulation [C]
   parameters%GDDS1     =     GDDS1_TABLE(CROPTYPE)    ! GDD from seeding to emergence
   parameters%GDDS2     =     GDDS2_TABLE(CROPTYPE)    ! GDD from seeding to initial vegetative 
   parameters%GDDS3     =     GDDS3_TABLE(CROPTYPE)    ! GDD from seeding to post vegetative 
   parameters%GDDS4     =     GDDS4_TABLE(CROPTYPE)    ! GDD from seeding to intial reproductive
   parameters%GDDS5     =     GDDS5_TABLE(CROPTYPE)    ! GDD from seeding to pysical maturity 
   parameters%C3C4      =      C3C4_TABLE(CROPTYPE)    ! photosynthetic pathway:  1. = c3 2. = c4
   parameters%AREF      =      AREF_TABLE(CROPTYPE)    ! reference maximum CO2 assimulation rate 
   parameters%PSNRF     =     PSNRF_TABLE(CROPTYPE)    ! CO2 assimulation reduction factor(0-1) (caused by non-modeling part,e.g.pest,weeds)
   parameters%I2PAR     =     I2PAR_TABLE(CROPTYPE)    ! Fraction of incoming solar radiation to photosynthetically active radiation
   parameters%TASSIM0   =   TASSIM0_TABLE(CROPTYPE)    ! Minimum temperature for CO2 assimulation [C]
   parameters%TASSIM1   =   TASSIM1_TABLE(CROPTYPE)    ! CO2 assimulation linearly increasing until temperature reaches T1 [C]
   parameters%TASSIM2   =   TASSIM2_TABLE(CROPTYPE)    ! CO2 assmilation rate remain at Aref until temperature reaches T2 [C]
   parameters%K         =         K_TABLE(CROPTYPE)    ! light extinction coefficient
   parameters%EPSI      =      EPSI_TABLE(CROPTYPE)    ! initial light use efficiency
   parameters%Q10MR     =     Q10MR_TABLE(CROPTYPE)    ! q10 for maintainance respiration
   parameters%FOLN_MX   =   FOLN_MX_TABLE(CROPTYPE)    ! foliage nitrogen concentration when f(n)=1 (%)
   parameters%LEFREEZ   =   LEFREEZ_TABLE(CROPTYPE)    ! characteristic T for leaf freezing [K]
   parameters%DILE_FC   =   DILE_FC_TABLE(CROPTYPE,:)  ! coeficient for temperature leaf stress death [1/s]
   parameters%DILE_FW   =   DILE_FW_TABLE(CROPTYPE,:)  ! coeficient for water leaf stress death [1/s]
   parameters%FRA_GR    =    FRA_GR_TABLE(CROPTYPE)    ! fraction of growth respiration
   parameters%LF_OVRC   =   LF_OVRC_TABLE(CROPTYPE,:)  ! fraction of leaf turnover  [1/s]
   parameters%ST_OVRC   =   ST_OVRC_TABLE(CROPTYPE,:)  ! fraction of stem turnover  [1/s]
   parameters%RT_OVRC   =   RT_OVRC_TABLE(CROPTYPE,:)  ! fraction of root tunrover  [1/s]
   parameters%LFMR25    =    LFMR25_TABLE(CROPTYPE)    ! leaf maintenance respiration at 25C [umol CO2/m**2  /s]
   parameters%STMR25    =    STMR25_TABLE(CROPTYPE)    ! stem maintenance respiration at 25C [umol CO2/kg bio/s]
   parameters%RTMR25    =    RTMR25_TABLE(CROPTYPE)    ! root maintenance respiration at 25C [umol CO2/kg bio/s]
   parameters%GRAINMR25 = GRAINMR25_TABLE(CROPTYPE)    ! grain maintenance respiration at 25C [umol CO2/kg bio/s]
   parameters%LFPT      =      LFPT_TABLE(CROPTYPE,:)  ! fraction of carbohydrate flux to leaf
   parameters%STPT      =      STPT_TABLE(CROPTYPE,:)  ! fraction of carbohydrate flux to stem
   parameters%RTPT      =      RTPT_TABLE(CROPTYPE,:)  ! fraction of carbohydrate flux to root
   parameters%GRAINPT   =   GRAINPT_TABLE(CROPTYPE,:)  ! fraction of carbohydrate flux to grain
   parameters%BIO2LAI   =   BIO2LAI_TABLE(CROPTYPE)    ! leaf are per living leaf biomass [m^2/kg]
  END IF

!------------------------------------------------------------------------------------------!
! Transfer global parameters
!------------------------------------------------------------------------------------------!

   parameters%CO2        =         CO2_TABLE
   parameters%O2         =          O2_TABLE
   parameters%TIMEAN     =      TIMEAN_TABLE
   parameters%FSATMX     =      FSATMX_TABLE
   parameters%Z0SNO      =       Z0SNO_TABLE
   parameters%SSI        =         SSI_TABLE
   parameters%SWEMX      =       SWEMX_TABLE
   parameters%RSURF_SNOW =  RSURF_SNOW_TABLE

! ----------------------------------------------------------------------
!  Transfer soil parameters
! ----------------------------------------------------------------------

    do isoil = 1, size(soiltype)
      parameters%BEXP(isoil)   = BEXP_TABLE   (SOILTYPE(isoil))
      parameters%DKSAT(isoil)  = DKSAT_TABLE  (SOILTYPE(isoil))
      parameters%DWSAT(isoil)  = DWSAT_TABLE  (SOILTYPE(isoil))
      parameters%PSISAT(isoil) = PSISAT_TABLE (SOILTYPE(isoil))
      parameters%QUARTZ(isoil) = QUARTZ_TABLE (SOILTYPE(isoil))
      parameters%SMCDRY(isoil) = SMCDRY_TABLE (SOILTYPE(isoil))
      parameters%SMCMAX(isoil) = SMCMAX_TABLE (SOILTYPE(isoil))
      parameters%SMCREF(isoil) = SMCREF_TABLE (SOILTYPE(isoil))
      parameters%SMCWLT(isoil) = SMCWLT_TABLE (SOILTYPE(isoil))
    end do
    
    parameters%F1     = F1_TABLE(SOILTYPE(1))
    parameters%REFDK  = REFDK_TABLE
    parameters%REFKDT = REFKDT_TABLE

! ----------------------------------------------------------------------
! Transfer GENPARM parameters
! ----------------------------------------------------------------------
    parameters%CSOIL  = CSOIL_TABLE
    parameters%ZBOT   = ZBOT_TABLE
    parameters%CZIL   = CZIL_TABLE

    FRZK   = FRZK_TABLE
    parameters%KDT    = parameters%REFKDT * parameters%DKSAT(1) / parameters%REFDK
    parameters%SLOPE  = SLOPE_TABLE(SLOPETYPE)

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
      parameters%FRZX = FRZK * FRZFACT
    END IF

    parameters%mxsnalb = 0.84
    parameters%mnsnalb = 0.55
    parameters%sndecayexp = 0.01

    parameters%t_ulimit = 2.5
    parameters%t_mlimit = 2.0
    parameters%t_llimit = 0.5
    parameters%snowf_scalef = 1.0    
 END SUBROUTINE TRANSFER_MP_PARAMETERS
