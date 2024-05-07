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
module noahmpglacier3911_Mod
!BOP
!
! !MODULE: noahmpglacier3911_Mod
!
! !DESCRIPTION:
!
! This module provides the definition of derived data type used to
! control the operation of noahmpglacier3911 model. It also provides the entry method
! for the initialization of noahmpglacier3911-specific variables. The derived
! data type {\tt noahmpglacier3911\_struc} includes the variables that specify
! the runtime options and other control variables as described below:
!
! \begin{description}
! \item[rfile]
!  name of the noahmpglacier3911 restart file
! \item[rformat]
!  format of restart file (binary or netcdf) for noahmpglacier3911
! \item[LDT\_ncvar\_vegetype]
!   LDT NetCDF variable name for land cover type index
! \item[LDT\_ncvar\_soiltype]
!   LDT NetCDF variable name for soil type index
! \item[LDT\_ncvar\_slopetype]
!   LDT NetCDF variable name for slope type for Noah baseflow
! \item[LDT\_ncvar\_tbot]
!   LDT NetCDF variable name for deep-layer soil temperature
! \item[LDT\_ncvar\_pblh]
!   LDT NetCDF variable name for planetary boundary layer height
! \item[ts]
!   noahmpglacier3911 model time step in second
! \item[count]
!   variable to keep track of the number of timesteps before an output
! \item[rstInterval]
!   restart writing interval
! \item[outInterval]
!   output writing interval
! \item[noahmp36]
!  noahmpglacier3911 model specific variables
! \item[forc\_count]
!   counter of forcing data
! \item[landuse\_tbl\_name]
!   Noah model landuse parameter table
! \item[soil\_tbl\_name]
!   Noah model soil parameter table
! \item[gen\_tbl\_name]
!   Noah model general parameter table
! \item[noahmp\_tbl\_name]
!   NoahMP parameter table
! \item[landuse\_scheme\_name]
!   landuse classification scheme
! \item[soil\_scheme\_name]
!   soil classification scheme
! \item[dveg\_opt]
!   vegetation model
! \item[crs\_opt]
!   canopy stomatal resistance
! \item[btr\_opt]
!   soil moisture factor for stomatal resistance
! \item[run\_opt]
!   runoff and groundwater
! \item[sfc\_opt]
!   surface layer drag coefficients (CH \& CM)
! \item[frz\_opt]
!   supercooled liquid water
! \item[inf\_opt]
!   frozen soil permeability
! \item[rad\_opt]
!   radiation transfer
! \item[alb\_opt]
!   snow surface albedo
! \item[snf\_opt]
!   rainfall \& snowfall
! \item[tbot\_opt]
!   lower boundary of soil temperature
! \item[stc\_opt]
!   snow/soil temperature time scheme
! \item[nslcats]
!   the number of total soil types in parameter table
! \item[nlucats]
!   the number of total land cover types in parameter table
! \item[nslpcats]
!   the number of total slope category for Noah baseflow
! \item[dt]
!   time step in seconds
! \item[nsoil]
!   number of soil layers
! \item[sldpth]
!   thickness of soil layers
! \item[nsnow]
!   maximum number of snow layers
! \item[urban\_vegetype]
!   urban land cover type index
! \item[ice\_flag]
!   ice flag: 0 = no ice, 1 = ice
! \item[st\_flag]
!   surface type 1=soil, 2=lake
! \item[sc\_idx]
!   soil color type
! \item[iz0tlnd]
!   option of Chen adjustment of Czil
! \item[zlvl]
!   reference height of temperature and humidity
! \end{description}
!
! !REVISION HISTORY:
!  This module is generated with the Model Implementation Toolkit developed
!  by Shugong Wang for the NASA Land Information System Version 7. The initial 
!  specification of the module is defined by Sujay Kumar. 
!  9/4/14: Shugong Wang Initial implementation for LIS 7 and noahmpglacier3911
!
! !USES:
    use noahmpglacier3911_module
    use LIS_constantsMod, only : LIS_CONST_PATH_LEN

    implicit none

    PRIVATE
    !-------------------------------------------------------------------------
    ! PUBLIC MEMBER FUNCTIONS
    !-------------------------------------------------------------------------
    public :: noahmpglacier3911_ini
    !-------------------------------------------------------------------------
    ! PUBLIC TYPES
    !-------------------------------------------------------------------------
    public :: noahmpgl3911_struc
!EOP
    type, public :: noahmpglacier_type_dec
       character(len=LIS_CONST_PATH_LEN) :: rfile
       character*256      :: rformat
       !-------------------------------------------------------------------------
       ! Parameter file names
       !-------------------------------------------------------------------------
       character*128      :: LDT_ncvar_shdfac_monthly
       character*128      :: LDT_ncvar_vegetype
       character*128      :: LDT_ncvar_soiltype
       character*128      :: LDT_ncvar_slopetype
       character*128      :: LDT_ncvar_smceq
       character*128      :: LDT_ncvar_tbot
       character*128      :: LDT_ncvar_pblh
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
       real               :: init_albold
       real               :: init_sneqvo
       real, pointer      :: init_stc(:)
       real, pointer      :: init_sh2o(:)
       real, pointer      :: init_smc(:)
       real               :: init_tah
       real               :: init_eah
       real               :: init_fwet
       real               :: init_canliq
       real               :: init_canice
       real               :: init_tv
       real               :: init_tg(1)
       real               :: init_qsnow
       !integer            :: init_isnow
       !real, pointer      :: init_zss(:)
       real               :: init_snowh(1)
       real               :: init_sneqv(1)
       !real, pointer      :: init_snowice(:)
       !real, pointer      :: init_snowliq(:)
       real               :: init_zwt
       real               :: init_wa
       real               :: init_wt
       real               :: init_wslake
       real               :: init_lfmass
       real               :: init_rtmass
       real               :: init_stmass
       real               :: init_wood
       real               :: init_stblcp
       real               :: init_fastcp
       real               :: init_lai
       real               :: init_sai
       real               :: init_cm
       real               :: init_ch
       real               :: init_tauss
       real               :: init_smcwtd
       real               :: init_deeprech
       real               :: init_rech
       real               :: init_zlvl 
       !-------------------------------------------------------------------------
       ! Constant Parameter
       !-------------------------------------------------------------------------
       character(len=256) :: landuse_tbl_name
       character(len=256) :: soil_tbl_name
       character(len=256) :: gen_tbl_name
       character(len=256) :: noahmp_tbl_name
       character(len=256) :: landuse_scheme_name
       character(len=256) :: soil_scheme_name
       integer            :: alb_opt
       integer            :: snf_opt
       integer            :: tbot_opt
       integer            :: stc_opt
       integer            :: gla_opt
       integer            :: nslcats
       integer            :: nlucats
       integer            :: nslpcats
       real               :: dt
       integer            :: nsoil
       real, pointer      :: sldpth(:)
       integer            :: nsnow
       integer            :: urban_vegetype
       integer            :: ice_flag
       integer            :: st_flag
       integer            :: sc_idx
       integer            :: iz0tlnd
       !real               :: zlvl
       type(noahmpgldec), pointer :: noahmpgl(:)
    end type noahmpglacier_type_dec
    
    type(noahmpglacier_type_dec), pointer :: noahmpgl3911_struc(:)
 
contains 

!BOP
!
! !ROUTINE: noahmpglacier3911_ini
! \label{noahmpglacier3911_ini}
!
! !INTERFACE:
  subroutine noahmpglacier3911_ini()
! !USES:
    use LIS_coreMod, only : LIS_rc
    use LIS_logMod, only : LIS_verify
    use LIS_timeMgrMod, only : LIS_clock,  LIS_calendar, &
         LIS_update_timestep, LIS_registerAlarm
    use LIS_surfaceModelDataMod, only : LIS_sfmodel_struc
! !DESCRIPTION:
!
!  This routine creates the datatypes and allocates memory for noahmpglacier3911-specific
!  variables. It also invokes the routine to read the runtime specific options
!  for noahmpglacier3911 from the configuration file.
!
!  The routines invoked are:
!  \begin{description}
!   \item[noahmpglacier3911\_readcrd](\ref{noahmpglacier3911_readcrd}) \newline
!    reads the runtime options for noahmpglacier3911 model
!  \end{description}
!EOP
    implicit none        
    integer  :: n, t     
    integer  :: status   

    ! allocate memory for nest 
    allocate(noahmpgl3911_struc(LIS_rc%nnest))
 
    ! read configuation information from lis.config file
    call noahmpglacier3911_readcrd()
    
    do n=1, LIS_rc%nnest
       ! allocate memory for all tiles in current nest 
       allocate(noahmpgl3911_struc(n)%noahmpgl(LIS_rc%npatch(n, LIS_rc%glacier_index)))
       !------------------------------------------------------------------------
       ! allocate memory for vector variables passed to model interfaces        
       ! TODO: check the following allocation statements carefully!
       !------------------------------------------------------------------------
       ! allocate memory for multilevel spatial parameter
#if 0
       do t=1, LIS_rc%npatch(n, LIS_rc%glacier_index)
          allocate(noahmpgl3911_struc(n)%noahmpgl(t)%shdfac_monthly(12))
          allocate(noahmpgl3911_struc(n)%noahmpgl(t)%smceq(noahmpgl3911_struc(n)%nsoil))
       enddo
#endif

! allocate memory for state variables
       do t=1, LIS_rc%npatch(n, LIS_rc%glacier_index)
          allocate(noahmpgl3911_struc(n)%noahmpgl(t)%sstc(&
               noahmpgl3911_struc(n)%nsoil + noahmpgl3911_struc(n)%nsnow))
          allocate(noahmpgl3911_struc(n)%noahmpgl(t)%sh2o(&
               noahmpgl3911_struc(n)%nsoil))
          allocate(noahmpgl3911_struc(n)%noahmpgl(t)%smc(&
               noahmpgl3911_struc(n)%nsoil))
          allocate(noahmpgl3911_struc(n)%noahmpgl(t)%zss(&
               noahmpgl3911_struc(n)%nsoil + noahmpgl3911_struc(n)%nsnow))
          allocate(noahmpgl3911_struc(n)%noahmpgl(t)%snowice(&
               noahmpgl3911_struc(n)%nsnow))
          allocate(noahmpgl3911_struc(n)%noahmpgl(t)%snowliq(&
               noahmpgl3911_struc(n)%nsnow))
       enddo
       ! initialize forcing variables to zeros
       do t=1, LIS_rc%npatch(n, LIS_rc%glacier_index)
          noahmpgl3911_struc(n)%noahmpgl(t)%lwdown = 0.0
          noahmpgl3911_struc(n)%noahmpgl(t)%swdown = 0.0
          noahmpgl3911_struc(n)%noahmpgl(t)%psurf = 0.0
          noahmpgl3911_struc(n)%noahmpgl(t)%prcp = 0.0
          noahmpgl3911_struc(n)%noahmpgl(t)%tair = 0.0
          noahmpgl3911_struc(n)%noahmpgl(t)%qair = 0.0
          noahmpgl3911_struc(n)%noahmpgl(t)%wind_e = 0.0
          noahmpgl3911_struc(n)%noahmpgl(t)%wind_n = 0.0
       enddo ! end of tile (t) loop
!------------------------------------------------------------------------
! Model timestep Alarm
!------------------------------------------------------------------------
       noahmpgl3911_struc(n)%forc_count = 0
       
       call LIS_update_timestep(LIS_rc, n, noahmpgl3911_struc(n)%ts)
       
       call LIS_registerAlarm("noahmpglacier3911 model alarm",&
            noahmpgl3911_struc(n)%ts, &
            noahmpgl3911_struc(n)%ts)
       
       call LIS_registerAlarm("noahmpglacier3911 restart alarm", &
            noahmpgl3911_struc(n)%ts,&
            noahmpgl3911_struc(n)%rstInterval)

#if 0 
!------------------------------------------------------------------------
! TODO: setup number of soil moisture/temperature layers and depth here  
!------------------------------------------------------------------------
       ! TODO: set number of soil moisture layers in surface model
       LIS_sfmodel_struc(n)%nsm_layers = noahmpgl3911_struc(n)%nsoil
       ! TODO: set number of soil temperature layers in surface model
       LIS_sfmodel_struc(n)%nst_layers = noahmpgl3911_struc(n)%nsoil
       allocate(LIS_sfmodel_struc(n)%lyrthk(noahmpgl3911_struc(n)%nsoil))
       LIS_sfmodel_struc(n)%lyrthk(:) = noahmpgl3911_struc(n)%sldpth(:) 
#endif
       LIS_sfmodel_struc(n)%ts = noahmpgl3911_struc(n)%ts
    enddo

  end subroutine noahmpglacier3911_ini
end module noahmpglacier3911_Mod
