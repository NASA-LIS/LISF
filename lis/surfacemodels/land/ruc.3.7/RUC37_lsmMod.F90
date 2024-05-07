!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!#include "LIS_misc.h"
module RUC37_lsmMod
!BOP
!
! !MODULE: RUC37_lsmMod
!
! !DESCRIPTION:
!
! This module provides the definition of derived data type used to
! control the operation of RUC37 model. It also provides the entry method
! for the initialization of RUC37-specific variables. The derived
! data type {\tt RUC37\_struc} includes the variables that specify
! the runtime options and other control variables as described below:
!
! \begin{description}
! \item[rfile]
!  name of the RUC37 restart file
! \item[rformat]
!  format of restart file (binary or netcdf) for RUC37
! \item[LDT\_ncvar\_vegetype]
!   LDT NetCDF variable name for vegetation category
! \item[LDT\_ncvar\_soiltype]
!   LDT NetCDF variable name for soil category
! \item[LDT\_ncvar\_albbck]
!   LDT NetCDF variable name for background snow-free albedo (0.0-1.0).
! \item[LDT\_ncvar\_tbot]
!   LDT NetCDF variable name for deep-soil time-invariant temperature (k).  representing sort of a mean annual air temperature.
! \item[LDT\_ncvar\_snoalb]
!   LDT NetCDF variable name for maximum snow albedo over deep snow (0.0-1.0)
! \item[ts]
!   RUC37 model time step in second
! \item[count]
!   variable to keep track of the number of timesteps before an output
! \item[rstInterval]
!   restart writing interval
! \item[outInterval]
!   output writing interval
! \item[ruc37]
!  RUC37 model specific variables
! \item[forc\_count]
!   counter of forcing data
! \item[dt]
!   time step (seconds).
! \item[soil\_layer\_thickness]
!   thicknesses of each soil level (m)
! \item[use\_local\_param]
!   .true. to use table values for albbck, shdfac, and z0brd; .false. to use values for albbck, shdfac, and z0brd as set in this driver routine
! \item[use\_2d\_lai\_map]
!   if rdlai2d == .true., then the xlai value that we pass to lsmruc will be used. if rdlai2d == .false., then xlai will be computed within lsmruc, from table minimum and maximum values in vegparm.tbl, and the current green vegetation fraction.
! \item[use\_monthly\_albedo\_map]
!   if usemonalb == .true., then the alb value passed to lsmruc will be used as the background snow-free albedo term.  if usemonalb == .false., then alb will be computed within lsmruc from minimum and maximum values in vegparm.tbl, and the current green vegetation fraction.
! \item[option\_iz0tlnd]
!   option to turn on (iz0tlnd=1) or off (iz0tlnd=0) the vegetation-category-dependent calculation of the zilitinkivich coefficient czil in the sfcdif subroutines.
! \item[option\_sfcdif]
!   option to use previous (sfcdif\_option=0) or updated (sfcdif\_option=1) version of sfcdif subroutine.
! \item[landuse\_tbl\_name]
!   noah model landuse parameter table
! \item[soil\_tbl\_name]
!   noah model soil parameter table
! \item[landuse\_scheme\_name]
!   landuse classification scheme
! \item[soil\_scheme\_name]
!   soil classification scheme
! \item[nsoil]
!   number of soil levels.
! \item[water\_class\_num]
!   number of water category in llanduse classification
! \item[ice\_class\_num]
!   number of ice category in llanduse classification
! \end{description}
!
! !REVISION HISTORY:
!  This module is generated with the Model Implementation Toolkit developed
!  by Shugong Wang for the NASA Land Information System Version 7. The initial 
!  specification of the module is defined by Sujay Kumar. 
!  1/15/15: Shugong Wang Initial implementation for LIS 7 and RUC37
!
! !USES:
    use RUC37_module
    use LIS_constantsMod, only : LIS_CONST_PATH_LEN

    implicit none

    PRIVATE
    !-------------------------------------------------------------------------
    ! PUBLIC MEMBER FUNCTIONS
    !-------------------------------------------------------------------------
    public :: RUC37_ini
    !-------------------------------------------------------------------------
    ! PUBLIC TYPES
    !-------------------------------------------------------------------------
    public :: RUC37_struc
!EOP
    type, public :: RUC37_type_dec
        character(len=LIS_CONST_PATH_LEN) :: rfile
        character*256      :: rformat
        !-------------------------------------------------------------------------
        ! Parameter file names
        !-------------------------------------------------------------------------
        character*128      :: LDT_ncvar_vegetype
        character*128      :: LDT_ncvar_soiltype
        character*128      :: LDT_ncvar_albedo_monthly
        character*128      :: LDT_ncvar_shdfac_monthly
        character*128      :: LDT_ncvar_z0brd_monthly
        character*128      :: LDT_ncvar_lai_monthly
        character*128      :: LDT_ncvar_albbck
        character*128      :: LDT_ncvar_tbot
        character*128      :: LDT_ncvar_snoalb
        !-------------------------------------------------------------------------
        ! ts, Count, rstInterval, outInterval
        !-------------------------------------------------------------------------
        real               :: ts
        integer            :: count
        real               :: rstInterval
        integer            :: outInterval
        integer            :: forc_count
        !-------------------------------------------------------------------------
        ! Initial Model State for cold start
        !-------------------------------------------------------------------------
        real               :: init_emiss
        real               :: init_ch
        real               :: init_cm
        real               :: init_sneqv
        real               :: init_snowh
        real               :: init_snowc
        real               :: init_canwat
        real               :: init_alb
        real, pointer      :: init_smc(:)
        real, pointer      :: init_sho(:)
        real, pointer      :: init_stc(:)
        real, pointer      :: init_smfr(:)
        real, pointer      :: init_keepfr(:)
        real               :: init_tskin
        real               :: init_qvg
        real               :: init_qsfc
        real               :: init_qcg
        real               :: init_qsg
        real               :: init_snt75cm
        real               :: init_tsnav
        real               :: init_soilm
        real               :: init_smroot
        !-------------------------------------------------------------------------
        ! Constant Parameter
        !-------------------------------------------------------------------------
        real               :: dt
        real, pointer      :: soil_layer_thickness(:)
        real               :: zlvl
        real               :: zlvl_wind 
        logical            :: use_local_param
        logical            :: use_2d_lai_map
        logical            :: use_monthly_albedo_map
        integer            :: option_iz0tlnd
        integer            :: option_sfcdif
        character(len=256) :: landuse_tbl_name
        character(len=256) :: soil_tbl_name
        character(len=256) :: gen_tbl_name
        character(len=256) :: landuse_scheme_name
        character(len=256) :: soil_scheme_name
        integer            :: nsoil
        integer            :: water_class_num
        integer            :: ice_class_num
        integer            :: urban_class_num
        type(RUC37dec), pointer :: ruc37(:)
    end type RUC37_type_dec

    type(RUC37_type_dec), pointer :: RUC37_struc(:)
 
contains 

!BOP
!
! !ROUTINE: RUC37_ini
! \label{RUC37_ini}
!
! !INTERFACE:
    subroutine RUC37_ini()
! !USES:
        use LIS_coreMod, only : LIS_rc
        use LIS_logMod, only : LIS_verify
        use LIS_timeMgrMod, only : LIS_clock,  LIS_calendar, &
            LIS_update_timestep, LIS_registerAlarm
        use LIS_surfaceModelDataMod, only : LIS_sfmodel_struc
! !DESCRIPTION:
!
!  This routine creates the datatypes and allocates memory for RUC37-specific
!  variables. It also invokes the routine to read the runtime specific options
!  for RUC37 from the configuration file.
!
!  The routines invoked are:
!  \begin{description}
!   \item[RUC37\_readcrd](\ref{RUC37_readcrd}) \newline
!    reads the runtime options for RUC37 model
!  \end{description}
!EOP
        implicit none        
        integer  :: n, t     
        integer  :: status   

        ! allocate memory for nest 
        allocate(RUC37_struc(LIS_rc%nnest))
 
        ! read configuation information from lis.config file
        call RUC37_readcrd()

        do n=1, LIS_rc%nnest
            ! allocate memory for all tiles in current nest 
            allocate(RUC37_struc(n)%ruc37(LIS_rc%npatch(n, LIS_rc%lsm_index)))
            !------------------------------------------------------------------------
            ! allocate memory for vector variables passed to model interfaces        
            ! TODO: check the following allocation statements carefully!
            !------------------------------------------------------------------------
            ! allocate memory for state variables
            do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
                allocate(RUC37_struc(n)%ruc37(t)%smc(RUC37_struc(n)%nsoil))
                allocate(RUC37_struc(n)%ruc37(t)%sho(RUC37_struc(n)%nsoil))
                allocate(RUC37_struc(n)%ruc37(t)%stc(RUC37_struc(n)%nsoil))
                allocate(RUC37_struc(n)%ruc37(t)%smfr(RUC37_struc(n)%nsoil))
                allocate(RUC37_struc(n)%ruc37(t)%keepfr(RUC37_struc(n)%nsoil))
                allocate(RUC37_struc(n)%ruc37(t)%albedo_monthly(12))
                allocate(RUC37_struc(n)%ruc37(t)%shdfac_monthly(12))
                allocate(RUC37_struc(n)%ruc37(t)%z0brd_monthly(12))
                allocate(RUC37_struc(n)%ruc37(t)%lai_monthly(12))
            enddo

            ! initialize forcing variables to zeros
            do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
                RUC37_struc(n)%ruc37(t)%lwdown = 0.0
                RUC37_struc(n)%ruc37(t)%swdown = 0.0
                RUC37_struc(n)%ruc37(t)%psurf = 0.0
                RUC37_struc(n)%ruc37(t)%rainf = 0.0
                RUC37_struc(n)%ruc37(t)%snowf = 0.0
                RUC37_struc(n)%ruc37(t)%tair = 0.0
                RUC37_struc(n)%ruc37(t)%qair = 0.0
                RUC37_struc(n)%ruc37(t)%wind_e = 0.0
                RUC37_struc(n)%ruc37(t)%wind_n = 0.0
            enddo ! end of tile (t) loop

            ! intialize run total variables to zero
            do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
                RUC37_struc(n)%ruc37(t)%snowfallac = 0.0
                RUC37_struc(n)%ruc37(t)%acsnow     = 0.0  
                RUC37_struc(n)%ruc37(t)%sfcevp     = 0.0 
            enddo
            !------------------------------------------------------------------------
            ! Model timestep Alarm
            !------------------------------------------------------------------------
            RUC37_struc(n)%forc_count = 0

            call LIS_update_timestep(LIS_rc, n, RUC37_struc(n)%ts)

            call LIS_registerAlarm("RUC37 model alarm",&
                                   RUC37_struc(n)%ts, &
                                   RUC37_struc(n)%ts)

            call LIS_registerAlarm("RUC37 restart alarm", &
                                   RUC37_struc(n)%ts,&
                                   RUC37_struc(n)%rstInterval)
            !------------------------------------------------------------------------
            ! TODO: setup number of soil moisture/temperature layers and depth here  
            !------------------------------------------------------------------------
            ! TODO: set number of soil moisture layers in surface model
            LIS_sfmodel_struc(n)%nsm_layers = RUC37_struc(n)%nsoil
            ! TODO: set number of soil temperature layers in surface model
            LIS_sfmodel_struc(n)%nst_layers = RUC37_struc(n)%nsoil
            allocate(LIS_sfmodel_struc(n)%lyrthk(RUC37_struc(n)%nsoil))
            LIS_sfmodel_struc(n)%ts = RUC37_struc(n)%ts
        enddo
    end subroutine RUC37_ini
end module RUC37_lsmMod
