!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LDT_misc.h"
module CLM45_parmsMod
!BOP
!
! !MODULE: CLM45_parmsMod
!
! !DESCRIPTION:
!  The code in this file implements routines
!  to read CLM-4.5 parameter data. 
!  \subsubsection{Overview}
!  This routines in this module provides routines
!  to read the CLM-4.5 parameter file data.
!
! !REVISION HISTORY:
!
!  08 Nov 2016: H. Beaudoing: Initial specification for CLM-4.5 model
!
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif

  use ESMF
  use LDT_coreMod
  use LDT_historyMod
  use LDT_paramDataMod
  use LDT_logMod
  use LDT_constantsMod, only : LDT_CONST_PATH_LEN
  use LDT_paramMaskCheckMod

  implicit none

  PRIVATE
!------------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!------------------------------------------------------------------------------
  public :: CLM45parms_init    !allocates memory for required structures
  public :: CLM45parms_writeHeader
  public :: CLM45parms_writeData

!------------------------------------------------------------------------------
! !PUBLIC TYPES:
!------------------------------------------------------------------------------
  public :: CLM45_struc

  type, public :: clm45_type_dec

     real      :: clm45_undef
     character(len=LDT_CONST_PATH_LEN) :: domainfile   ! CLM-4.5 precomputed domain file
     character(len=LDT_CONST_PATH_LEN) :: surfacefile  ! CLM-4.5 precomputed surface file
     character*50  :: clm45parms_gridtransform
     character*50  :: clm45parms_proj
     real          :: clm45parms_gridDesc(20)
     character(len=32)    :: clm45_param_mode

! -  CLM-specific:
     integer       :: ncorners  
     integer       :: nlevsoi      ! number of soil layers
     integer       :: lsmpft       ! number of pft 
     integer       :: numurbl      ! number of urban types
     integer       :: nlevurb      ! urban levels
     integer       :: numrad       ! number of radiation
     integer       :: nglcec       ! number of glacier layers
     integer       :: nglcecp1     ! glacier melt layers
     type(LDT_paramEntry) :: clm45
! - domain.vars.txt variables -
     type(LDT_paramEntry) :: xc    !        double xc(nj, ni) ;
     type(LDT_paramEntry) :: yc    !        double yc(nj, ni) ;
     type(LDT_paramEntry) :: xv    !        double xv(nj, ni, nv) ;
     type(LDT_paramEntry) :: yv    !        double yv(nj, ni, nv) ;
     type(LDT_paramEntry) :: mask  !        int mask(nj, ni) ;
     type(LDT_paramEntry) :: area  !        double area(nj, ni) ;
     type(LDT_paramEntry) :: frac  !        double frac(nj, ni) ;
!! - surfdata.vars.txt variables -
     type(LDT_paramEntry) :: mxsoil_color
     type(LDT_paramEntry) :: SOIL_COLOR !   int SOIL_COLOR(lsmlat, lsmlon) ;
     type(LDT_paramEntry) :: PCT_SAND !double PCT_SAND(nlevsoi, lsmlat, lsmlon);
     type(LDT_paramEntry) :: PCT_CLAY !double PCT_CLAY(nlevsoi, lsmlat, lsmlon);
     type(LDT_paramEntry) :: ORGANIC  !double ORGANIC(nlevsoi, lsmlat, lsmlon);
     type(LDT_paramEntry) :: FMAX     !        double FMAX(lsmlat, lsmlon) ;
     type(LDT_paramEntry) :: LANDFRAC_PFT !double LANDFRAC_PFT(lsmlat, lsmlon) ;
     type(LDT_paramEntry) :: PFTDATA_MASK !int PFTDATA_MASK(lsmlat, lsmlon) ;
     type(LDT_paramEntry) :: PCT_PFT  !double PCT_PFT(lsmpft, lsmlat, lsmlon) ;
     type(LDT_paramEntry) :: MONTHLY_LAI !double MONTHLY_LAI(time, lsmpft, lsmlat, lsmlon) ;
     type(LDT_paramEntry) :: MONTHLY_SAI !double MONTHLY_SAI(time, lsmpft, lsmlat, lsmlon) ;
     type(LDT_paramEntry) :: MONTHLY_HEIGHT_TOP !double MONTHLY_HEIGHT_TOP(time, lsmpft, lsmlat, lsmlon) ;
     type(LDT_paramEntry) :: MONTHLY_HEIGHT_BOT !double MONTHLY_HEIGHT_BOT(time, lsmpft, lsmlat, lsmlon) ;
     type(LDT_paramEntry) :: time !        int time(time) ;
     type(LDT_paramEntry) :: SRFAREA !        double AREA(lsmlat, lsmlon) ;
     type(LDT_paramEntry) :: LONGXY !        double LONGXY(lsmlat, lsmlon) ;
     type(LDT_paramEntry) :: LATIXY !        double LATIXY(lsmlat, lsmlon) ;
     type(LDT_paramEntry) :: EF1_BTR !        double EF1_BTR(lsmlat, lsmlon) ;
     type(LDT_paramEntry) :: EF1_FET !        double EF1_FET(lsmlat, lsmlon) ;
     type(LDT_paramEntry) :: EF1_FDT !        double EF1_FDT(lsmlat, lsmlon) ;
     type(LDT_paramEntry) :: EF1_SHR !        double EF1_SHR(lsmlat, lsmlon) ;
     type(LDT_paramEntry) :: EF1_GRS !        double EF1_GRS(lsmlat, lsmlon) ;
     type(LDT_paramEntry) :: EF1_CRP !        double EF1_CRP(lsmlat, lsmlon) ;
     type(LDT_paramEntry) :: CANYON_HWR !        double CANYON_HWR(numurbl, lsmlat, lsmlon) ;
     type(LDT_paramEntry) :: EM_IMPROAD !        double EM_IMPROAD(numurbl, lsmlat, lsmlon) ;
     type(LDT_paramEntry) :: EM_PERROAD !        double EM_PERROAD(numurbl, lsmlat, lsmlon) ;
     type(LDT_paramEntry) :: EM_ROOF  !double EM_ROOF(numurbl, lsmlat, lsmlon) ;
     type(LDT_paramEntry) :: EM_WALL  !double EM_WALL(numurbl, lsmlat, lsmlon) ;
     type(LDT_paramEntry) :: HT_ROOF  !double HT_ROOF(numurbl, lsmlat, lsmlon) ;
     type(LDT_paramEntry) :: THICK_ROOF !double THICK_ROOF(numurbl, lsmlat, lsmlon) ;
     type(LDT_paramEntry) :: THICK_WALL !double THICK_WALL(numurbl, lsmlat, lsmlon) ;
     type(LDT_paramEntry) :: T_BUILDING_MAX !        double T_BUILDING_MAX(numurbl, lsmlat, lsmlon) ;
     type(LDT_paramEntry) :: T_BUILDING_MIN !        double T_BUILDING_MIN(numurbl, lsmlat, lsmlon) ;
     type(LDT_paramEntry) :: WIND_HGT_CANYON !       double WIND_HGT_CANYON(numurbl, lsmlat, lsmlon) ;
     type(LDT_paramEntry) :: WTLUNIT_ROOF    !       double WTLUNIT_ROOF(numurbl, lsmlat, lsmlon) ;
     type(LDT_paramEntry) :: WTROAD_PERV !        double WTROAD_PERV(numurbl, lsmlat, lsmlon) ;
     type(LDT_paramEntry) :: ALB_IMPROAD_DIR !        double ALB_IMPROAD_DIR(numrad, numurbl, lsmlat, lsmlon) ;
     type(LDT_paramEntry) :: ALB_IMPROAD_DIF !        double ALB_IMPROAD_DIF(numrad, numurbl, lsmlat, lsmlon) ;
     type(LDT_paramEntry) :: ALB_PERROAD_DIR !        double ALB_PERROAD_DIR(numrad, numurbl, lsmlat, lsmlon) ;
     type(LDT_paramEntry) :: ALB_PERROAD_DIF !        double ALB_PERROAD_DIF(numrad, numurbl, lsmlat, lsmlon) ;
     type(LDT_paramEntry) :: ALB_ROOF_DIR    !        double ALB_ROOF_DIR(numrad, numurbl, lsmlat, lsmlon) ;
     type(LDT_paramEntry) :: ALB_ROOF_DIF    !        double ALB_ROOF_DIF(numrad, numurbl, lsmlat, lsmlon) ;
     type(LDT_paramEntry) :: ALB_WALL_DIR    !        double ALB_WALL_DIR(numrad, numurbl, lsmlat, lsmlon) ;
     type(LDT_paramEntry) :: ALB_WALL_DIF    !        double ALB_WALL_DIF(numrad, numurbl, lsmlat, lsmlon) ;
     type(LDT_paramEntry) :: TK_ROOF !        double TK_ROOF(nlevurb, numurbl, lsmlat, lsmlon) ;
     type(LDT_paramEntry) :: TK_WALL !        double TK_WALL(nlevurb, numurbl, lsmlat, lsmlon) ;
     type(LDT_paramEntry) :: TK_IMPROAD !        double TK_IMPROAD(nlevurb, numurbl, lsmlat, lsmlon) ;
     type(LDT_paramEntry) :: CV_ROOF !        double CV_ROOF(nlevurb, numurbl, lsmlat, lsmlon) ;
     type(LDT_paramEntry) :: CV_WALL !        double CV_WALL(nlevurb, numurbl, lsmlat, lsmlon) ;
     type(LDT_paramEntry) :: CV_IMPROAD !        double CV_IMPROAD(nlevurb, numurbl, lsmlat, lsmlon) ;
     type(LDT_paramEntry) :: NLEV_IMPROAD !        int NLEV_IMPROAD(numurbl, lsmlat, lsmlon) ;
     type(LDT_paramEntry) :: peatf !        double peatf(lsmlat, lsmlon) ;
     type(LDT_paramEntry) :: abm   !        int abm(lsmlat, lsmlon) ;
     type(LDT_paramEntry) :: gdp   !        double gdp(lsmlat, lsmlon) ;
     type(LDT_paramEntry) :: SLOPE !        double SLOPE(lsmlat, lsmlon) ;
     type(LDT_paramEntry) :: STD_ELEV !        double STD_ELEV(lsmlat, lsmlon) ;
     type(LDT_paramEntry) :: binfl !        double binfl(lsmlat, lsmlon) ;
     type(LDT_paramEntry) :: Ws    !        double Ws(lsmlat, lsmlon) ;
     type(LDT_paramEntry) :: Dsmax !        double Dsmax(lsmlat, lsmlon) ;
     type(LDT_paramEntry) :: Ds    !        double Ds(lsmlat, lsmlon) ;
     type(LDT_paramEntry) :: LAKEDEPTH !    double LAKEDEPTH(lsmlat, lsmlon) ;
     type(LDT_paramEntry) :: F0    !        double F0(lsmlat, lsmlon) ;
     type(LDT_paramEntry) :: P3    !        double P3(lsmlat, lsmlon) ;
     type(LDT_paramEntry) :: ZWT0  !        double ZWT0(lsmlat, lsmlon) ;
     type(LDT_paramEntry) :: PCT_WETLAND !  double PCT_WETLAND(lsmlat, lsmlon) ;
     type(LDT_paramEntry) :: PCT_LAKE    !  double PCT_LAKE(lsmlat, lsmlon) ;
     type(LDT_paramEntry) :: PCT_GLACIER !  double PCT_GLACIER(lsmlat, lsmlon) ;
     type(LDT_paramEntry) :: GLC_MEC     !  double GLC_MEC(nglcecp1) ;
     type(LDT_paramEntry) :: PCT_GLC_MEC !  double PCT_GLC_MEC(nglcec, lsmlat, lsmlon) ;
     type(LDT_paramEntry) :: PCT_GLC_MEC_GIC      !        double PCT_GLC_MEC_GIC(nglcec, lsmlat, lsmlon) ;
     type(LDT_paramEntry) :: PCT_GLC_MEC_ICESHEET !        double PCT_GLC_MEC_ICESHEET(nglcec, lsmlat, lsmlon) ;
     type(LDT_paramEntry) :: PCT_GLC_GIC !  double PCT_GLC_GIC(lsmlat, lsmlon) ;
     type(LDT_paramEntry) :: PCT_GLC_ICESHEET !        double PCT_GLC_ICESHEET(lsmlat, lsmlon) ;
     type(LDT_paramEntry) :: TOPO_GLC_MEC     !        double TOPO_GLC_MEC(nglcec, lsmlat, lsmlon) ;
     type(LDT_paramEntry) :: TOPO  !        double TOPO(lsmlat, lsmlon) ;
     type(LDT_paramEntry) :: PCT_URBAN !double PCT_URBAN(numurbl, lsmlat, lsmlon) ;
     type(LDT_paramEntry) :: URBAN_REGION_ID !        int URBAN_REGION_ID(lsmlat, lsmlon) ;
  end type clm45_type_dec

  type(clm45_type_dec), allocatable :: CLM45_struc(:)

contains

!BOP
! 
! !ROUTINE: CLM45parms_init
! \label{CLM45parms_init}
! 
! !INTERFACE:
  subroutine CLM45parms_init(flag)

! !USES:
   use LDT_logMod,  only : LDT_verify, LDT_endrun, &
                           LDT_getNextUnitNumber, LDT_releaseUnitNumber
   use LDT_fileIOMod,only : LDT_readDomainConfigSpecs
   use LDT_paramOptCheckMod, only: LDT_gridOptChecks
!
! !DESCRIPTION:
!
! Community Land Model (CLM) v4.5 model parameters.
! 
!  The routines invoked are: 
!  \begin{description}
!   \item[clm45parmssetup](\ref{clm45parmssetup}) \newline
!    calls the registry to invoke the clm45parms setup methods. 
!  \end{description}
!
!EOP
   implicit none
   integer  :: flag
   integer  :: n
   integer  :: c,r,m,k
   integer  :: rc
   integer  :: file_status
   logical  :: file_exists
   logical  :: clm45_select 
   real, allocatable          :: clm45parms_gridDesc(:,:)
   character*50, allocatable  :: clm45parms_gridtransform(:)
   character*50               :: clm45parms_proj
   character(len=32)          :: clm45_param_mode

 ! _____________________________________________________________________

   allocate(CLM45_struc(LDT_rc%nnest))
   clm45_select = .false.
   do n=1,LDT_rc%nnest
      ! - CLM-4.5 parameters:
      call set_param_attribs(CLM45_struc(n)%clm45,"CLM45","","-")
      if( CLM45_struc(n)%clm45%selectOpt == 1 ) then
         clm45_select = .true.
         LDT_rc%monthlyData(n) = .true.  ! turn on 12-mo time for output
      endif
   enddo

   if( clm45_select ) then
      
    write(LDT_logunit,*)" - - - - - - - - - - CLM-4.5 LSM Parameters - - - - - - - - - - - - -"
    write(LDT_logunit,*)" ** MSG:  Only CLM-4.5 Model Parameters available."

    allocate(clm45parms_gridDesc(LDT_rc%nnest,20))
    allocate(clm45parms_gridtransform(LDT_rc%nnest))

    do n=1,LDT_rc%nnest
     CLM45_struc(n)%ncorners = 4
     CLM45_struc(n)%nlevsoi = 10
     CLM45_struc(n)%lsmpft = 17
     CLM45_struc(n)%numurbl = 3
     CLM45_struc(n)%numrad = 2
     CLM45_struc(n)%nlevurb = 5
     CLM45_struc(n)%nglcecp1 = 11
     CLM45_struc(n)%nglcec = 10

     call set_param_attribs(CLM45_struc(n)%xc,"xc", &
          long_name="longitude of grid cell center",units="degrees_east")
     call set_param_attribs(CLM45_struc(n)%yc,"yc", &
          long_name="latitude of grid cell center",units="degrees_north")
     call set_param_attribs(CLM45_struc(n)%xv,"xv", &
          vlevels=CLM45_struc(n)%ncorners, &
          long_name="longitude of grid cell verticies",units="degrees_east")
     call set_param_attribs(CLM45_struc(n)%yv,"yv", &
          vlevels=CLM45_struc(n)%ncorners, &
          long_name="latitude of grid cell verticies",units="degrees_north")
! LANDMASK is already defined -> set as LDT_LSMparam_struc(n)%landmask?
!     call set_param_attribs(CLM45_struc(n)%mask,"LANDMASK", &
     call set_param_attribs(CLM45_struc(n)%mask,"CLMMASK", &
          long_name="domain mask",units="unitless")
     call set_param_attribs(CLM45_struc(n)%area,"area", &
          long_name="area of grid cell in radians squared",units="radian2")
     call set_param_attribs(CLM45_struc(n)%frac,"frac", &
          long_name="fraction of grid cell that is active",units="unitless")

     call set_param_attribs(CLM45_struc(n)%mxsoil_color,"mxsoil_color", &
          long_name="maximum numbers of soil colors",units="unitless")
     call set_param_attribs(CLM45_struc(n)%SOIL_COLOR,"SOIL_COLOR", &
          long_name="soil color",units="unitless")
     call set_param_attribs(CLM45_struc(n)%PCT_SAND,"PCT_SAND", &
          vlevels=CLM45_struc(n)%nlevsoi, &
          long_name="percent sand",units="unitless")
     call set_param_attribs(CLM45_struc(n)%PCT_CLAY,"PCT_CLAY", &
          vlevels=CLM45_struc(n)%nlevsoi, &
          long_name="percent clay",units="unitless")
     call set_param_attribs(CLM45_struc(n)%ORGANIC,"ORGANIC", &
          vlevels=CLM45_struc(n)%nlevsoi, &
          long_name="organic matter density at soil levels",  &
          units="kg/m3 (assumed carbon content 0.58 gC per gOM)")
     call set_param_attribs(CLM45_struc(n)%FMAX,"FMAX", &
          long_name="maximum fractional saturated area",units="unitless")
     call set_param_attribs(CLM45_struc(n)%LANDFRAC_PFT,"LANDFRAC_PFT", &
          long_name="land fraction from pft dataset",units="unitless")
     call set_param_attribs(CLM45_struc(n)%PFTDATA_MASK,"PFTDATA_MASK", &
          long_name="land mask from pft dataset, indicative of real/fake points",units="unitless")
     call set_param_attribs(CLM45_struc(n)%PCT_PFT,"PCT_PFT", &
          vlevels=CLM45_struc(n)%lsmpft, &
          long_name="percent plant functional type of gridcell",units="unitless")
     call set_param_attribs(CLM45_struc(n)%MONTHLY_LAI,"MONTHLY_LAI", &
          vlevels=CLM45_struc(n)%lsmpft, &
          zlevels=12, &
          long_name="monthly leaf area index",units="unitless")
     call set_param_attribs(CLM45_struc(n)%MONTHLY_SAI,"MONTHLY_SAI", &
          vlevels=CLM45_struc(n)%lsmpft, &
          zlevels=12, &
          long_name="monthly stem area index",units="unitless")
     call set_param_attribs(CLM45_struc(n)%MONTHLY_HEIGHT_TOP,"MONTHLY_HEIGHT_TOP", &
          vlevels=CLM45_struc(n)%lsmpft, &
          zlevels=12, &
          long_name="monthly height top",units="meters")
     call set_param_attribs(CLM45_struc(n)%MONTHLY_HEIGHT_BOT,"MONTHLY_HEIGHT_BOT", &
          vlevels=CLM45_struc(n)%lsmpft, &
          zlevels=12, &
          long_name="monthly height bottom",units="meters")
     call set_param_attribs(CLM45_struc(n)%time,"time", &
          long_name="Calendar month",units="month")
     call set_param_attribs(CLM45_struc(n)%SRFAREA,"AREA", &
          long_name="area",units="km^2")
     call set_param_attribs(CLM45_struc(n)%LONGXY,"LONGXY", &
          long_name="longitude",units="degrees east")
     call set_param_attribs(CLM45_struc(n)%LATIXY,"LATIXY", &
          long_name="latitude",units="degrees north")
     call set_param_attribs(CLM45_struc(n)%EF1_BTR,"EF1_BTR", &
          long_name="EF btr (isoprene)",units="unitless")
     call set_param_attribs(CLM45_struc(n)%EF1_FET,"EF1_FET", &
          long_name="EF fet (isoprene)",units="unitless")
     call set_param_attribs(CLM45_struc(n)%EF1_FDT,"EF1_FDT", &
          long_name="EF fdt (isoprene)",units="unitless")
     call set_param_attribs(CLM45_struc(n)%EF1_SHR,"EF1_SHR", &
          long_name="EF shr (isoprene)",units="unitless")
     call set_param_attribs(CLM45_struc(n)%EF1_GRS,"EF1_GRS", &
          long_name="EF grs (isoprene)",units="unitless")
     call set_param_attribs(CLM45_struc(n)%EF1_CRP,"EF1_CRP", &
          long_name="EF crp (isoprene)",units="unitless")
     call set_param_attribs(CLM45_struc(n)%CANYON_HWR,"CANYON_HWR", &
          vlevels=CLM45_struc(n)%numurbl, &
          long_name="canyon height to width ratio",units="unitless")
     call set_param_attribs(CLM45_struc(n)%EM_IMPROAD,"EM_IMPROAD", &
          vlevels=CLM45_struc(n)%numurbl, &
          long_name="emissivity of impervious road",units="unitless")
     call set_param_attribs(CLM45_struc(n)%EM_PERROAD,"EM_PERROAD", &
          vlevels=CLM45_struc(n)%numurbl, &
          long_name="emissivity of pervious road",units="unitless")
     call set_param_attribs(CLM45_struc(n)%EM_ROOF,"EM_ROOF", &
          vlevels=CLM45_struc(n)%numurbl, &
          long_name="emissivity of roof",units="unitless")
     call set_param_attribs(CLM45_struc(n)%EM_WALL,"EM_WALL", &
          vlevels=CLM45_struc(n)%numurbl, &
          long_name="emissivity of wall",units="unitless")
     call set_param_attribs(CLM45_struc(n)%HT_ROOF,"HT_ROOF", &
          vlevels=CLM45_struc(n)%numurbl, &
          long_name="height of roof",units="unitless")
     call set_param_attribs(CLM45_struc(n)%THICK_ROOF,"THICK_ROOF", &
          vlevels=CLM45_struc(n)%numurbl, &
          long_name="thickness of roof",units="meters")
     call set_param_attribs(CLM45_struc(n)%THICK_WALL,"THICK_WALL", &
          vlevels=CLM45_struc(n)%numurbl, &
          long_name="thickness of wall",units="meters")
     call set_param_attribs(CLM45_struc(n)%T_BUILDING_MAX,"T_BUILDING_MAX", &
          vlevels=CLM45_struc(n)%numurbl, &
          long_name="maximum interior building temperature",units="K")
     call set_param_attribs(CLM45_struc(n)%T_BUILDING_MIN,"T_BUILDING_MIN", &
          vlevels=CLM45_struc(n)%numurbl, &
          long_name="miniimum interior building temperature",units="K")
     call set_param_attribs(CLM45_struc(n)%WIND_HGT_CANYON,"WIND_HGT_CANYON", &
          vlevels=CLM45_struc(n)%numurbl, &
          long_name="height of wind in canyon",units="meters")
     call set_param_attribs(CLM45_struc(n)%WTLUNIT_ROOF,"WTLUNIT_ROOF", &
          vlevels=CLM45_struc(n)%numurbl, &
          long_name="fraction of roof",units="unitless")
     call set_param_attribs(CLM45_struc(n)%WTROAD_PERV,"WTROAD_PERV", &
          vlevels=CLM45_struc(n)%numurbl, &
          long_name="fraction of pervious road",units="unitless")
     call set_param_attribs(CLM45_struc(n)%ALB_IMPROAD_DIR,"ALB_IMPROAD_DIR", &
          vlevels=CLM45_struc(n)%numurbl, &
          zlevels=CLM45_struc(n)%numrad, &
          long_name="direct albedo of impervious road",units="unitless")
     call set_param_attribs(CLM45_struc(n)%ALB_IMPROAD_DIF,"ALB_IMPROAD_DIF", &
          vlevels=CLM45_struc(n)%numurbl, &
          zlevels=CLM45_struc(n)%numrad, &
          long_name="diffuse albedo of impervious road",units="unitless")
     call set_param_attribs(CLM45_struc(n)%ALB_PERROAD_DIR,"ALB_PERROAD_DIR", &
          vlevels=CLM45_struc(n)%numurbl, &
          zlevels=CLM45_struc(n)%numrad, &
          long_name="direct albedo of pervious road",units="unitless")
     call set_param_attribs(CLM45_struc(n)%ALB_PERROAD_DIF,"ALB_PERROAD_DIF", &
          vlevels=CLM45_struc(n)%numurbl, &
          zlevels=CLM45_struc(n)%numrad, &
          long_name="diffuse albedo of pervious road",units="unitless")
     call set_param_attribs(CLM45_struc(n)%ALB_ROOF_DIR,"ALB_ROOF_DIR", &
          vlevels=CLM45_struc(n)%numurbl, &
          zlevels=CLM45_struc(n)%numrad, &
          long_name="direct albedo of roof",units="unitless")
     call set_param_attribs(CLM45_struc(n)%ALB_ROOF_DIF,"ALB_ROOF_DIF", &
          vlevels=CLM45_struc(n)%numurbl, &
          zlevels=CLM45_struc(n)%numrad, &
          long_name="diffuse albedo of roof",units="unitless")
     call set_param_attribs(CLM45_struc(n)%ALB_WALL_DIR,"ALB_WALL_DIR", &
          vlevels=CLM45_struc(n)%numurbl, &
          zlevels=CLM45_struc(n)%numrad, &
          long_name="direct albedo of wall",units="unitless")
     call set_param_attribs(CLM45_struc(n)%ALB_WALL_DIF,"ALB_WALL_DIF", &
          vlevels=CLM45_struc(n)%numurbl, &
          zlevels=CLM45_struc(n)%numrad, &
          long_name="diffuse albedo of wall",units="unitless")
     call set_param_attribs(CLM45_struc(n)%TK_ROOF,"TK_ROOF", &
          vlevels=CLM45_struc(n)%numurbl, &
          zlevels=CLM45_struc(n)%nlevurb, &
          long_name="thermal conductivity of roof",units="W/m*K")
     call set_param_attribs(CLM45_struc(n)%TK_WALL,"TK_WALL", &
          vlevels=CLM45_struc(n)%numurbl, &
          zlevels=CLM45_struc(n)%nlevurb, &
          long_name="thermal conductivity of wall",units="W/m*K")
     call set_param_attribs(CLM45_struc(n)%TK_IMPROAD,"TK_IMPROAD", &
          vlevels=CLM45_struc(n)%numurbl, &
          zlevels=CLM45_struc(n)%nlevurb, &
          long_name="thermal conductivity of impervious road",units="W/m*K")
     call set_param_attribs(CLM45_struc(n)%CV_ROOF,"CV_ROOF", &
          vlevels=CLM45_struc(n)%numurbl, &
          zlevels=CLM45_struc(n)%nlevurb, &
          long_name="volumetric heat capacity of roof",units="J/m^3*K")
     call set_param_attribs(CLM45_struc(n)%CV_WALL,"CV_WALL", &
          vlevels=CLM45_struc(n)%numurbl, &
          zlevels=CLM45_struc(n)%nlevurb, &
          long_name="volumetric heat capacity of wall",units="J/m^3*K")
     call set_param_attribs(CLM45_struc(n)%CV_IMPROAD,"CV_IMPROAD", &
          vlevels=CLM45_struc(n)%numurbl, &
          zlevels=CLM45_struc(n)%nlevurb, &
          long_name="volumetric heat capacity of impervious road",units="J/m^3*K")
     call set_param_attribs(CLM45_struc(n)%NLEV_IMPROAD,"NLEV_IMPROAD", &
          vlevels=CLM45_struc(n)%numurbl, &
          long_name="number of impervious road layers",units="unitless")
     call set_param_attribs(CLM45_struc(n)%peatf,"peatf", &
          long_name="peatland fraction",units="unitless")
     call set_param_attribs(CLM45_struc(n)%abm,"abm", &
          long_name="agricultural fire peak month",units="unitless")
     call set_param_attribs(CLM45_struc(n)%gdp,"gdp", &
          long_name="gdp",units="unitless")
     call set_param_attribs(CLM45_struc(n)%SLOPE,"SLOPE", &
          long_name="mean topographic slope",units="degrees")
     call set_param_attribs(CLM45_struc(n)%STD_ELEV,"STD_ELEV", &
          long_name="standard deviation of elevation",units="m")
     call set_param_attribs(CLM45_struc(n)%binfl,"binfl", &
          long_name="VIC b parameter for the Variable Infiltration Capacity Curve",units="unitless")
     call set_param_attribs(CLM45_struc(n)%Ws,"Ws", &
          long_name="VIC Ws parameter for the ARNO Curve",units="unitless")
     call set_param_attribs(CLM45_struc(n)%Dsmax,"Dsmax", &
          long_name="VIC Dsmax parameter for the ARNO curve",units="mm/day")
     call set_param_attribs(CLM45_struc(n)%Ds,"Ds", &
          long_name="VIC Ds parameter for the ARNO curve",units="unitless")
     call set_param_attribs(CLM45_struc(n)%LAKEDEPTH,"LAKEDEPTH", &
          long_name="lake depth",units="m")
     call set_param_attribs(CLM45_struc(n)%F0,"F0", &
          long_name="maximum gridcell fractional inundated area",units="unitless")
     call set_param_attribs(CLM45_struc(n)%P3,"P3", &
          long_name="coefficient for qflx_surf_lag for finundated",units="s/mm")
     call set_param_attribs(CLM45_struc(n)%ZWT0,"ZWT0", &
          long_name="decay factor for finundated",units="m")
     call set_param_attribs(CLM45_struc(n)%PCT_WETLAND,"PCT_WETLAND", &
          long_name="percent wetland",units="unitless")
     call set_param_attribs(CLM45_struc(n)%PCT_LAKE,"PCT_LAKE", &
          long_name="percent lake",units="unitless")
     call set_param_attribs(CLM45_struc(n)%PCT_GLACIER,"PCT_GLACIER", &
          long_name="percent glacier",units="unitless")
     call set_param_attribs(CLM45_struc(n)%GLC_MEC,"GLC_MEC", &
          long_name="Glacier elevation class",units="m")
     call set_param_attribs(CLM45_struc(n)%PCT_GLC_MEC,"PCT_GLC_MEC", &
          vlevels=CLM45_struc(n)%nglcec, &
          long_name="percent glacier for each glacier elevation class",units="unitless")
     call set_param_attribs(CLM45_struc(n)%PCT_GLC_MEC_GIC,"PCT_GLC_MEC_GIC", &
          vlevels=CLM45_struc(n)%nglcec, &
          long_name="percent smaller glaciers and ice caps for each glacier elevation class",units="unitless")
     call set_param_attribs(CLM45_struc(n)%PCT_GLC_MEC_ICESHEET,"PCT_GLC_MEC_ICESHEET", &
          vlevels=CLM45_struc(n)%nglcec, &
          long_name="percent ice sheet for each glacier elevation class",units="unitless")
     call set_param_attribs(CLM45_struc(n)%PCT_GLC_GIC,"PCT_GLC_GIC", &
          long_name="percent ice caps/glaciers",units="unitless")
     call set_param_attribs(CLM45_struc(n)%PCT_GLC_ICESHEET,"PCT_GLC_ICESHEET", &
          long_name="percent ice sheet",units="unitless")
     call set_param_attribs(CLM45_struc(n)%TOPO_GLC_MEC,"TOPO_GLC_MEC", &
          vlevels=CLM45_struc(n)%nglcec, &
          long_name="mean elevation on glacier elevation classes",units="m")
     call set_param_attribs(CLM45_struc(n)%TOPO,"TOPO", &
          long_name="mean elevation on land",units="m")
     call set_param_attribs(CLM45_struc(n)%PCT_URBAN,"PCT_URBAN", &
          vlevels=CLM45_struc(n)%numurbl, &
          long_name="percent urban for each density type",units="unitless")
     call set_param_attribs(CLM45_struc(n)%URBAN_REGION_ID, &
          "URBAN_REGION_ID", &
          long_name="urban region ID",units="unitless")

     allocate(CLM45_struc(n)%xc%dvalue(LDT_rc%lnc(n),LDT_rc%lnr(n),1,1))
     allocate(CLM45_struc(n)%yc%dvalue(LDT_rc%lnc(n),LDT_rc%lnr(n),1,1))
     allocate(CLM45_struc(n)%xv%dvalue(LDT_rc%lnc(n),LDT_rc%lnr(n),&
              CLM45_struc(n)%ncorners,1))
     allocate(CLM45_struc(n)%yv%dvalue(LDT_rc%lnc(n),LDT_rc%lnr(n),&
              CLM45_struc(n)%ncorners,1))
     allocate(CLM45_struc(n)%mask%dvalue(LDT_rc%lnc(n),LDT_rc%lnr(n),1,1))
     allocate(CLM45_struc(n)%area%dvalue(LDT_rc%lnc(n),LDT_rc%lnr(n),1,1))
     allocate(CLM45_struc(n)%frac%dvalue(LDT_rc%lnc(n),LDT_rc%lnr(n),1,1))

     allocate(CLM45_struc(n)%mxsoil_color%dvalue(1,1,1,1))
     allocate(CLM45_struc(n)%SOIL_COLOR%dvalue(LDT_rc%lnc(n),LDT_rc%lnr(n),1,1))
     allocate(CLM45_struc(n)%PCT_SAND%dvalue(LDT_rc%lnc(n),LDT_rc%lnr(n),CLM45_struc(n)%nlevsoi,1))
     allocate(CLM45_struc(n)%PCT_CLAY%dvalue(LDT_rc%lnc(n),LDT_rc%lnr(n),CLM45_struc(n)%nlevsoi,1))
     allocate(CLM45_struc(n)%ORGANIC%dvalue(LDT_rc%lnc(n),LDT_rc%lnr(n),CLM45_struc(n)%nlevsoi,1))
     allocate(CLM45_struc(n)%FMAX%dvalue(LDT_rc%lnc(n),LDT_rc%lnr(n),1,1))
     allocate(CLM45_struc(n)%LANDFRAC_PFT%dvalue(LDT_rc%lnc(n),LDT_rc%lnr(n),1,1))
     allocate(CLM45_struc(n)%PFTDATA_MASK%dvalue(LDT_rc%lnc(n),LDT_rc%lnr(n),1,1))
     allocate(CLM45_struc(n)%PCT_PFT%dvalue(LDT_rc%lnc(n),LDT_rc%lnr(n),CLM45_struc(n)%lsmpft,1))
     allocate(CLM45_struc(n)%MONTHLY_LAI%dvalue(LDT_rc%lnc(n),LDT_rc%lnr(n),CLM45_struc(n)%lsmpft,12))
     allocate(CLM45_struc(n)%MONTHLY_SAI%dvalue(LDT_rc%lnc(n),LDT_rc%lnr(n),CLM45_struc(n)%lsmpft,12))
     allocate(CLM45_struc(n)%MONTHLY_HEIGHT_TOP%dvalue(LDT_rc%lnc(n),LDT_rc%lnr(n),CLM45_struc(n)%lsmpft,12))
     allocate(CLM45_struc(n)%MONTHLY_HEIGHT_BOT%dvalue(LDT_rc%lnc(n),LDT_rc%lnr(n),CLM45_struc(n)%lsmpft,12))
     allocate(CLM45_struc(n)%time%dvalue(12,1,1,1))
     allocate(CLM45_struc(n)%SRFAREA%dvalue(LDT_rc%lnc(n),LDT_rc%lnr(n),1,1))
     allocate(CLM45_struc(n)%LONGXY%dvalue(LDT_rc%lnc(n),LDT_rc%lnr(n),1,1))
     allocate(CLM45_struc(n)%LATIXY%dvalue(LDT_rc%lnc(n),LDT_rc%lnr(n),1,1))
     allocate(CLM45_struc(n)%EF1_BTR%dvalue(LDT_rc%lnc(n),LDT_rc%lnr(n),1,1))
     allocate(CLM45_struc(n)%EF1_FET%dvalue(LDT_rc%lnc(n),LDT_rc%lnr(n),1,1))
     allocate(CLM45_struc(n)%EF1_FDT%dvalue(LDT_rc%lnc(n),LDT_rc%lnr(n),1,1))
     allocate(CLM45_struc(n)%EF1_SHR%dvalue(LDT_rc%lnc(n),LDT_rc%lnr(n),1,1))
     allocate(CLM45_struc(n)%EF1_GRS%dvalue(LDT_rc%lnc(n),LDT_rc%lnr(n),1,1))
     allocate(CLM45_struc(n)%EF1_CRP%dvalue(LDT_rc%lnc(n),LDT_rc%lnr(n),1,1))
     allocate(CLM45_struc(n)%CANYON_HWR%dvalue(LDT_rc%lnc(n),LDT_rc%lnr(n),CLM45_struc(n)%numurbl,1))
     allocate(CLM45_struc(n)%EM_IMPROAD%dvalue(LDT_rc%lnc(n),LDT_rc%lnr(n),CLM45_struc(n)%numurbl,1))
     allocate(CLM45_struc(n)%EM_PERROAD%dvalue(LDT_rc%lnc(n),LDT_rc%lnr(n),CLM45_struc(n)%numurbl,1))
     allocate(CLM45_struc(n)%EM_ROOF%dvalue(LDT_rc%lnc(n),LDT_rc%lnr(n),CLM45_struc(n)%numurbl,1))
     allocate(CLM45_struc(n)%EM_WALL%dvalue(LDT_rc%lnc(n),LDT_rc%lnr(n),CLM45_struc(n)%numurbl,1))
     allocate(CLM45_struc(n)%HT_ROOF%dvalue(LDT_rc%lnc(n),LDT_rc%lnr(n),CLM45_struc(n)%numurbl,1))
     allocate(CLM45_struc(n)%THICK_ROOF%dvalue(LDT_rc%lnc(n),LDT_rc%lnr(n),CLM45_struc(n)%numurbl,1))
     allocate(CLM45_struc(n)%THICK_WALL%dvalue(LDT_rc%lnc(n),LDT_rc%lnr(n),CLM45_struc(n)%numurbl,1))
     allocate(CLM45_struc(n)%T_BUILDING_MAX%dvalue(LDT_rc%lnc(n),LDT_rc%lnr(n),CLM45_struc(n)%numurbl,1))
     allocate(CLM45_struc(n)%T_BUILDING_MIN%dvalue(LDT_rc%lnc(n),LDT_rc%lnr(n),CLM45_struc(n)%numurbl,1))
     allocate(CLM45_struc(n)%WIND_HGT_CANYON%dvalue(LDT_rc%lnc(n),LDT_rc%lnr(n),CLM45_struc(n)%numurbl,1))
     allocate(CLM45_struc(n)%WTLUNIT_ROOF%dvalue(LDT_rc%lnc(n),LDT_rc%lnr(n),CLM45_struc(n)%numurbl,1))
     allocate(CLM45_struc(n)%WTROAD_PERV%dvalue(LDT_rc%lnc(n),LDT_rc%lnr(n),CLM45_struc(n)%numurbl,1))
     allocate(CLM45_struc(n)%ALB_IMPROAD_DIR%dvalue(LDT_rc%lnc(n),LDT_rc%lnr(n),CLM45_struc(n)%numurbl,CLM45_struc(n)%numrad))
     allocate(CLM45_struc(n)%ALB_IMPROAD_DIF%dvalue(LDT_rc%lnc(n),LDT_rc%lnr(n),CLM45_struc(n)%numurbl,CLM45_struc(n)%numrad))
     allocate(CLM45_struc(n)%ALB_PERROAD_DIR%dvalue(LDT_rc%lnc(n),LDT_rc%lnr(n),CLM45_struc(n)%numurbl,CLM45_struc(n)%numrad))
     allocate(CLM45_struc(n)%ALB_PERROAD_DIF%dvalue(LDT_rc%lnc(n),LDT_rc%lnr(n),CLM45_struc(n)%numurbl,CLM45_struc(n)%numrad))
     allocate(CLM45_struc(n)%ALB_ROOF_DIR%dvalue(LDT_rc%lnc(n),LDT_rc%lnr(n),CLM45_struc(n)%numurbl,CLM45_struc(n)%numrad))
     allocate(CLM45_struc(n)%ALB_ROOF_DIF%dvalue(LDT_rc%lnc(n),LDT_rc%lnr(n),CLM45_struc(n)%numurbl,CLM45_struc(n)%numrad))
     allocate(CLM45_struc(n)%ALB_WALL_DIR%dvalue(LDT_rc%lnc(n),LDT_rc%lnr(n),CLM45_struc(n)%numurbl,CLM45_struc(n)%numrad))
     allocate(CLM45_struc(n)%ALB_WALL_DIF%dvalue(LDT_rc%lnc(n),LDT_rc%lnr(n),CLM45_struc(n)%numurbl,CLM45_struc(n)%numrad))
     allocate(CLM45_struc(n)%TK_ROOF%dvalue(LDT_rc%lnc(n),LDT_rc%lnr(n),CLM45_struc(n)%numurbl,CLM45_struc(n)%nlevurb))
     allocate(CLM45_struc(n)%TK_WALL%dvalue(LDT_rc%lnc(n),LDT_rc%lnr(n),CLM45_struc(n)%numurbl,CLM45_struc(n)%nlevurb))
     allocate(CLM45_struc(n)%TK_IMPROAD%dvalue(LDT_rc%lnc(n),LDT_rc%lnr(n),CLM45_struc(n)%numurbl,CLM45_struc(n)%nlevurb))
     allocate(CLM45_struc(n)%CV_ROOF%dvalue(LDT_rc%lnc(n),LDT_rc%lnr(n),CLM45_struc(n)%numurbl,CLM45_struc(n)%nlevurb))
     allocate(CLM45_struc(n)%CV_WALL%dvalue(LDT_rc%lnc(n),LDT_rc%lnr(n),CLM45_struc(n)%numurbl,CLM45_struc(n)%nlevurb))
     allocate(CLM45_struc(n)%CV_IMPROAD%dvalue(LDT_rc%lnc(n),LDT_rc%lnr(n),CLM45_struc(n)%numurbl,CLM45_struc(n)%nlevurb))
     allocate(CLM45_struc(n)%NLEV_IMPROAD%dvalue(LDT_rc%lnc(n),LDT_rc%lnr(n),CLM45_struc(n)%numurbl,1))
     allocate(CLM45_struc(n)%peatf%dvalue(LDT_rc%lnc(n),LDT_rc%lnr(n),1,1))
     allocate(CLM45_struc(n)%abm%dvalue(LDT_rc%lnc(n),LDT_rc%lnr(n),1,1))
     allocate(CLM45_struc(n)%gdp%dvalue(LDT_rc%lnc(n),LDT_rc%lnr(n),1,1))
     allocate(CLM45_struc(n)%SLOPE%dvalue(LDT_rc%lnc(n),LDT_rc%lnr(n),1,1))
     allocate(CLM45_struc(n)%STD_ELEV%dvalue(LDT_rc%lnc(n),LDT_rc%lnr(n),1,1))
     allocate(CLM45_struc(n)%binfl%dvalue(LDT_rc%lnc(n),LDT_rc%lnr(n),1,1))
     allocate(CLM45_struc(n)%Ws%dvalue(LDT_rc%lnc(n),LDT_rc%lnr(n),1,1))
     allocate(CLM45_struc(n)%Dsmax%dvalue(LDT_rc%lnc(n),LDT_rc%lnr(n),1,1))
     allocate(CLM45_struc(n)%Ds%dvalue(LDT_rc%lnc(n),LDT_rc%lnr(n),1,1))
     allocate(CLM45_struc(n)%LAKEDEPTH%dvalue(LDT_rc%lnc(n),LDT_rc%lnr(n),1,1))
     allocate(CLM45_struc(n)%F0%dvalue(LDT_rc%lnc(n),LDT_rc%lnr(n),1,1))
     allocate(CLM45_struc(n)%P3%dvalue(LDT_rc%lnc(n),LDT_rc%lnr(n),1,1))
     allocate(CLM45_struc(n)%ZWT0%dvalue(LDT_rc%lnc(n),LDT_rc%lnr(n),1,1))
     allocate(CLM45_struc(n)%PCT_WETLAND%dvalue(LDT_rc%lnc(n),LDT_rc%lnr(n),1,1))
     allocate(CLM45_struc(n)%PCT_LAKE%dvalue(LDT_rc%lnc(n),LDT_rc%lnr(n),1,1))
     allocate(CLM45_struc(n)%PCT_GLACIER%dvalue(LDT_rc%lnc(n),LDT_rc%lnr(n),1,1))
     allocate(CLM45_struc(n)%GLC_MEC%dvalue(CLM45_struc(n)%nglcecp1,1,1,1))
     allocate(CLM45_struc(n)%PCT_GLC_MEC%dvalue(LDT_rc%lnc(n),LDT_rc%lnr(n),CLM45_struc(n)%nglcec,1))
     allocate(CLM45_struc(n)%PCT_GLC_MEC_GIC%dvalue(LDT_rc%lnc(n),LDT_rc%lnr(n),CLM45_struc(n)%nglcec,1))
     allocate(CLM45_struc(n)%PCT_GLC_MEC_ICESHEET%dvalue(LDT_rc%lnc(n),LDT_rc%lnr(n),CLM45_struc(n)%nglcec,1))
     allocate(CLM45_struc(n)%PCT_GLC_GIC%dvalue(LDT_rc%lnc(n),LDT_rc%lnr(n),1,1))
     allocate(CLM45_struc(n)%PCT_GLC_ICESHEET%dvalue(LDT_rc%lnc(n),LDT_rc%lnr(n),1,1))
     allocate(CLM45_struc(n)%TOPO_GLC_MEC%dvalue(LDT_rc%lnc(n),LDT_rc%lnr(n),CLM45_struc(n)%nglcec,1))
     allocate(CLM45_struc(n)%TOPO%dvalue(LDT_rc%lnc(n),LDT_rc%lnr(n),1,1))
     allocate(CLM45_struc(n)%PCT_URBAN%dvalue(LDT_rc%lnc(n),LDT_rc%lnr(n),CLM45_struc(n)%numurbl,1))
     allocate(CLM45_struc(n)%URBAN_REGION_ID%dvalue(LDT_rc%lnc(n),LDT_rc%lnr(n),1,1))

    enddo   !n
   endif  ! clm45_select

!-- CLM-4.5 Parameter Config Entries: 
!-- mode: "readin" or "create"
    call ESMF_ConfigFindLabel(LDT_config,"CLM45 parameter mode:",rc=rc)
    call ESMF_ConfigGetAttribute(LDT_config,clm45_param_mode,rc=rc)
    call LDT_verify(rc,'CLM45 parameter mode: not specified')
!    print*,'clm45_param_mode:',clm45_param_mode
!-- spatial transform: "none" when "readin"; "average"or"neighbor" when "create"
    call ESMF_ConfigFindLabel(LDT_config, &
                              "CLM45 param spatial transform:",rc=rc)
    do n=1,LDT_rc%nnest
       call ESMF_ConfigGetAttribute(LDT_config,clm45parms_gridtransform(n),rc=rc)
       call LDT_verify(rc,'CLM45 param spatial transform: not specified')
       CLM45_struc(n)%clm45parms_gridtransform = clm45parms_gridtransform(n)
    enddo
    call ESMF_ConfigFindLabel(LDT_config,"CLM45 param map projection:",rc=rc)
    do n=1,LDT_rc%nnest
       call ESMF_ConfigGetAttribute(LDT_config,clm45parms_proj,rc=rc)
       call LDT_verify(rc,'CLM45 param map projection: not specified')
       CLM45_struc(n)%clm45parms_proj = clm45parms_proj
    enddo
    do n=1,LDT_rc%nnest
     call LDT_readDomainConfigSpecs("CLM45",clm45parms_proj,clm45parms_gridDesc)
!     print*,'readDomainConfigSpecs:',clm45parms_gridDesc
     if ( clm45parms_proj == "latlon" ) then
       call LDT_gridOptChecks( n, "CLM45", &
            clm45parms_gridtransform(n), &
            clm45parms_proj, &
            clm45parms_gridDesc(n,9) )
       CLM45_struc(n)%clm45parms_gridDesc = clm45parms_gridDesc(n,:)
     endif
    enddo

    if (trim(clm45_param_mode) .eq. "readin") then
     call ESMF_ConfigFindLabel(LDT_config,"CLM45 domain file:",rc=rc)
     do n=1,LDT_rc%nnest
        call ESMF_ConfigGetAttribute(LDT_config,&
            CLM45_struc(n)%domainfile,rc=rc)
        call LDT_verify(rc,'CLM45 domain file: not specified')
     enddo
     call ESMF_ConfigFindLabel(LDT_config,"CLM45 surface file:",rc=rc)
     do n=1,LDT_rc%nnest
        call ESMF_ConfigGetAttribute(LDT_config,&
             CLM45_struc(n)%surfacefile,rc=rc)
        call LDT_verify(rc,'CLM45 surface file: not specified')
     enddo
    else  ! "create"
        write(LDT_logunit,*) "[ERR] CLM-4.5 'create' option not supported yet"
        write(LDT_logunit,*) " Program stopping ..."
        call LDT_endrun
    endif   !clm45_param_mode

    if (trim(clm45_param_mode) .eq. "readin") then
     do n=1,LDT_rc%nnest
      call clm45_read_params(n)
      write(LDT_logunit,*) "** MSG: Reset Landcover and Surfacetype for CLM-4.5"
      call getclm45lcst(n)
     enddo
    else  ! "create"
     !call clm45_create_params
    endif

    if ( clm45_select ) then
     deallocate(clm45parms_gridDesc)
     deallocate(clm45parms_gridtransform)
    endif
 end subroutine CLM45parms_init

 subroutine clm45_read_params(n)
  use LDT_logMod
  use LDT_paramDataMod
  use LDT_coreMod
 
  implicit none
  integer :: n

!-- Read in CLM-4.5 Parameter Datasets:

      write(LDT_logunit,*) "Reading domain values from "//trim(CLM45_struc(n)%domainfile)
      call read_netcdf_4d(n,trim(CLM45_struc(n)%domainfile),"xc",2, &
          "nj","ni","none","none",LDT_rc%lnc(n),LDT_rc%lnr(n),1,1,&
          CLM45_struc(n)%xc%dvalue(:,:,1,1))     
      write(LDT_logunit,*) "Done reading xc values."
      call read_netcdf_4d(n,trim(CLM45_struc(n)%domainfile),"yc",2, &
          "nj","ni","none","none",LDT_rc%lnc(n),LDT_rc%lnr(n),1,1,&
          CLM45_struc(n)%yc%dvalue(:,:,1,1))     
      write(LDT_logunit,*) "Done reading yc values."
      call read_netcdf_4d_xvyv(n,trim(CLM45_struc(n)%domainfile),"xv",3, &
          "nj","ni","nv","none", &
          LDT_rc%lnc(n),LDT_rc%lnr(n),CLM45_struc(n)%ncorners,1,&
          CLM45_struc(n)%xv%dvalue(:,:,:,1))     
      write(LDT_logunit,*) "Done reading xv values."
      call read_netcdf_4d_xvyv(n,trim(CLM45_struc(n)%domainfile),"yv",3, &
          "nj","ni","nv","none", &
          LDT_rc%lnc(n),LDT_rc%lnr(n),CLM45_struc(n)%ncorners,1,&
          CLM45_struc(n)%yv%dvalue(:,:,:,1))     
      write(LDT_logunit,*) "Done reading yv values."
      call read_netcdf_4d(n,trim(CLM45_struc(n)%domainfile),"mask",2, &
          "nj","ni","none","none",LDT_rc%lnc(n),LDT_rc%lnr(n),1,1,&
          CLM45_struc(n)%mask%dvalue(:,:,1,1))     
      write(LDT_logunit,*) "Done reading mask values."
      call read_netcdf_4d(n,trim(CLM45_struc(n)%domainfile),"area",2, &
          "nj","ni","none","none",LDT_rc%lnc(n),LDT_rc%lnr(n),1,1,&
          CLM45_struc(n)%area%dvalue(:,:,1,1))     
      write(LDT_logunit,*) "Done reading area values."
      call read_netcdf_4d(n,trim(CLM45_struc(n)%domainfile),"frac",2, &
          "nj","ni","none","none",LDT_rc%lnc(n),LDT_rc%lnr(n),1,1,&
          CLM45_struc(n)%frac%dvalue(:,:,1,1))     
      write(LDT_logunit,*) "Done reading frac values."

      write(LDT_logunit,*) "Reading surface values from "//trim(CLM45_struc(n)%surfacefile)
      call read_netcdf_4d(n,trim(CLM45_struc(n)%surfacefile),"mxsoil_color",0, &
          "none ","none","none","none",1,1,1,1,&
          CLM45_struc(n)%mxsoil_color%dvalue(1,1,1,1))     
      write(LDT_logunit,*) "Done reading mxsoil_color values."
      call read_netcdf_4d(n,trim(CLM45_struc(n)%surfacefile),"SOIL_COLOR",2, &
          "lsmlat","lsmlon","none","none",LDT_rc%lnc(n),LDT_rc%lnr(n),1,1,&
          CLM45_struc(n)%SOIL_COLOR%dvalue(:,:,1,1))     
      write(LDT_logunit,*) "Done reading SOIL_COLOR values."
      call read_netcdf_4d(n,trim(CLM45_struc(n)%surfacefile),"PCT_SAND",3, &
          "nlevsoi","lsmlat","lsmlon","none", &
          LDT_rc%lnc(n),LDT_rc%lnr(n),CLM45_struc(n)%nlevsoi,1, &
          CLM45_struc(n)%PCT_SAND%dvalue(:,:,:,1))     
      write(LDT_logunit,*) "Done reading PCT_SAND values."
      call read_netcdf_4d(n,trim(CLM45_struc(n)%surfacefile),"PCT_CLAY",3, &
          "nlevsoi","lsmlat","lsmlon","none", &
          LDT_rc%lnc(n),LDT_rc%lnr(n),CLM45_struc(n)%nlevsoi,1, &
          CLM45_struc(n)%PCT_CLAY%dvalue(:,:,:,1))     
      write(LDT_logunit,*) "Done reading PCT_CLAY values."
      call read_netcdf_4d(n,trim(CLM45_struc(n)%surfacefile),"ORGANIC",3, &
          "nlevsoi","lsmlat","lsmlon","none", &
          LDT_rc%lnc(n),LDT_rc%lnr(n),CLM45_struc(n)%nlevsoi,1, &
          CLM45_struc(n)%ORGANIC%dvalue(:,:,:,1))     
      write(LDT_logunit,*) "Done reading ORGANIC values."
      call read_netcdf_4d(n,trim(CLM45_struc(n)%surfacefile),"FMAX",2, &
          "lsmlat","lsmlon","none","none",LDT_rc%lnc(n),LDT_rc%lnr(n),1,1,&
          CLM45_struc(n)%FMAX%dvalue(:,:,1,1))     
      write(LDT_logunit,*) "Done reading FMAX values."
      call read_netcdf_4d(n,trim(CLM45_struc(n)%surfacefile),"LANDFRAC_PFT",2, &
          "lsmlat","lsmlon","none","none",LDT_rc%lnc(n),LDT_rc%lnr(n),1,1,&
          CLM45_struc(n)%LANDFRAC_PFT%dvalue(:,:,1,1))     
      write(LDT_logunit,*) "Done reading LANDFRAC_PFT values."
      call read_netcdf_4d(n,trim(CLM45_struc(n)%surfacefile),"PFTDATA_MASK",2, &
          "lsmlat","lsmlon","none","none",LDT_rc%lnc(n),LDT_rc%lnr(n),1,1,&
          CLM45_struc(n)%PFTDATA_MASK%dvalue(:,:,1,1))     
      write(LDT_logunit,*) "Done reading PFTDATA_MASK values."
      call read_netcdf_4d(n,trim(CLM45_struc(n)%surfacefile),"PCT_PFT",3, &
          "lsmpft","lsmlat","lsmlon","none", &
          LDT_rc%lnc(n),LDT_rc%lnr(n),CLM45_struc(n)%lsmpft,1, &
          CLM45_struc(n)%PCT_PFT%dvalue(:,:,:,1))     
      write(LDT_logunit,*) "Done reading PCT_PFT values."
      call read_netcdf_4d(n,trim(CLM45_struc(n)%surfacefile),"MONTHLY_LAI",4, &
          "time","lsmpft","lsmlat","lsmlon", &
          LDT_rc%lnc(n),LDT_rc%lnr(n),CLM45_struc(n)%lsmpft,12, &
          CLM45_struc(n)%MONTHLY_LAI%dvalue(:,:,:,:))     
      write(LDT_logunit,*) "Done reading MONTHLY_LAI values."
      call read_netcdf_4d(n,trim(CLM45_struc(n)%surfacefile),"MONTHLY_SAI",4, &
          "time","lsmpft","lsmlat","lsmlon", &
          LDT_rc%lnc(n),LDT_rc%lnr(n),CLM45_struc(n)%lsmpft,12, &
          CLM45_struc(n)%MONTHLY_SAI%dvalue(:,:,:,:))     
      write(LDT_logunit,*) "Done reading MONTHLY_SAI values."
      call read_netcdf_4d(n,trim(CLM45_struc(n)%surfacefile),"MONTHLY_HEIGHT_TOP",4, &
          "time","lsmpft","lsmlat","lsmlon", &
          LDT_rc%lnc(n),LDT_rc%lnr(n),CLM45_struc(n)%lsmpft,12, &
          CLM45_struc(n)%MONTHLY_HEIGHT_TOP%dvalue(:,:,:,:))     
      write(LDT_logunit,*) "Done reading MONTHLY_HEIGHT_TOP values."
      call read_netcdf_4d(n,trim(CLM45_struc(n)%surfacefile),"MONTHLY_HEIGHT_BOT",4, &
          "time","lsmpft","lsmlat","lsmlon", &
          LDT_rc%lnc(n),LDT_rc%lnr(n),CLM45_struc(n)%lsmpft,12, &
          CLM45_struc(n)%MONTHLY_HEIGHT_BOT%dvalue(:,:,:,:))     
      write(LDT_logunit,*) "Done reading MONTHLY_HEIGHT_BOT values."
      call read_netcdf_4d(n,trim(CLM45_struc(n)%surfacefile),"time",1, &
          "time","none","none","none",12,1,1,1,&
          CLM45_struc(n)%time%dvalue(:,1,1,1))     
      write(LDT_logunit,*) "Done reading time values."
      call read_netcdf_4d(n,trim(CLM45_struc(n)%surfacefile),"AREA",2, &
          "lsmlat","lsmlon","none","none",LDT_rc%lnc(n),LDT_rc%lnr(n),1,1,&
          CLM45_struc(n)%SRFAREA%dvalue(:,:,1,1))     
      write(LDT_logunit,*) "Done reading AREA values."
      call read_netcdf_4d(n,trim(CLM45_struc(n)%surfacefile),"LONGXY",2, &
          "lsmlat","lsmlon","none","none",LDT_rc%lnc(n),LDT_rc%lnr(n),1,1,&
          CLM45_struc(n)%LONGXY%dvalue(:,:,1,1))     
      write(LDT_logunit,*) "Done reading LONGXY values."
      call read_netcdf_4d(n,trim(CLM45_struc(n)%surfacefile),"LATIXY",2, &
          "lsmlat","lsmlon","none","none",LDT_rc%lnc(n),LDT_rc%lnr(n),1,1,&
          CLM45_struc(n)%LATIXY%dvalue(:,:,1,1))     
      write(LDT_logunit,*) "Done reading LATIXY values."
      call read_netcdf_4d(n,trim(CLM45_struc(n)%surfacefile),"EF1_BTR",2, &
          "lsmlat","lsmlon","none","none",LDT_rc%lnc(n),LDT_rc%lnr(n),1,1,&
          CLM45_struc(n)%EF1_BTR%dvalue(:,:,1,1))     
      write(LDT_logunit,*) "Done reading EF1_BTR values."
      call read_netcdf_4d(n,trim(CLM45_struc(n)%surfacefile),"EF1_FET",2, &
          "lsmlat","lsmlon","none","none",LDT_rc%lnc(n),LDT_rc%lnr(n),1,1,&
          CLM45_struc(n)%EF1_FET%dvalue(:,:,1,1))     
      write(LDT_logunit,*) "Done reading EF1_FET values."
      call read_netcdf_4d(n,trim(CLM45_struc(n)%surfacefile),"EF1_FDT",2, &
          "lsmlat","lsmlon","none","none",LDT_rc%lnc(n),LDT_rc%lnr(n),1,1,&
          CLM45_struc(n)%EF1_FDT%dvalue(:,:,1,1))     
      write(LDT_logunit,*) "Done reading EF1_FDT values."
      call read_netcdf_4d(n,trim(CLM45_struc(n)%surfacefile),"EF1_SHR",2, &
          "lsmlat","lsmlon","none","none",LDT_rc%lnc(n),LDT_rc%lnr(n),1,1,&
          CLM45_struc(n)%EF1_SHR%dvalue(:,:,1,1))     
      write(LDT_logunit,*) "Done reading EF1_SHR values."
      call read_netcdf_4d(n,trim(CLM45_struc(n)%surfacefile),"EF1_GRS",2, &
          "lsmlat","lsmlon","none","none",LDT_rc%lnc(n),LDT_rc%lnr(n),1,1,&
          CLM45_struc(n)%EF1_GRS%dvalue(:,:,1,1))     
      write(LDT_logunit,*) "Done reading EF1_GRS values."
      call read_netcdf_4d(n,trim(CLM45_struc(n)%surfacefile),"EF1_CRP",2, &
          "lsmlat","lsmlon","none","none",LDT_rc%lnc(n),LDT_rc%lnr(n),1,1,&
          CLM45_struc(n)%EF1_CRP%dvalue(:,:,1,1))     
      write(LDT_logunit,*) "Done reading EF1_CRP values."
      call read_netcdf_4d(n,trim(CLM45_struc(n)%surfacefile),"CANYON_HWR",3, &
          "numurbl","lsmlat","lsmlon","none", &
          LDT_rc%lnc(n),LDT_rc%lnr(n),CLM45_struc(n)%numurbl,1, &
          CLM45_struc(n)%CANYON_HWR%dvalue(:,:,:,1))     
      write(LDT_logunit,*) "Done reading CANYON_HWR values."
      call read_netcdf_4d(n,trim(CLM45_struc(n)%surfacefile),"EM_IMPROAD",3, &
          "numurbl","lsmlat","lsmlon","none", &
          LDT_rc%lnc(n),LDT_rc%lnr(n),CLM45_struc(n)%numurbl,1, &
          CLM45_struc(n)%EM_IMPROAD%dvalue(:,:,:,1))     
      write(LDT_logunit,*) "Done reading EM_IMPROAD values."
      call read_netcdf_4d(n,trim(CLM45_struc(n)%surfacefile),"EM_PERROAD",3, &
          "numurbl","lsmlat","lsmlon","none", &
          LDT_rc%lnc(n),LDT_rc%lnr(n),CLM45_struc(n)%numurbl,1, &
          CLM45_struc(n)%EM_PERROAD%dvalue(:,:,:,1))     
      write(LDT_logunit,*) "Done reading EM_PERROAD values."
      call read_netcdf_4d(n,trim(CLM45_struc(n)%surfacefile),"EM_ROOF",3, &
          "numurbl","lsmlat","lsmlon","none", &
          LDT_rc%lnc(n),LDT_rc%lnr(n),CLM45_struc(n)%numurbl,1, &
          CLM45_struc(n)%EM_ROOF%dvalue(:,:,:,1))     
      write(LDT_logunit,*) "Done reading EM_ROOF values."
      call read_netcdf_4d(n,trim(CLM45_struc(n)%surfacefile),"EM_WALL",3, &
          "numurbl","lsmlat","lsmlon","none", &
          LDT_rc%lnc(n),LDT_rc%lnr(n),CLM45_struc(n)%numurbl,1, &
          CLM45_struc(n)%EM_WALL%dvalue(:,:,:,1))     
      write(LDT_logunit,*) "Done reading EM_WALL values."
      call read_netcdf_4d(n,trim(CLM45_struc(n)%surfacefile),"HT_ROOF",3, &
          "numurbl","lsmlat","lsmlon","none", &
          LDT_rc%lnc(n),LDT_rc%lnr(n),CLM45_struc(n)%numurbl,1, &
          CLM45_struc(n)%HT_ROOF%dvalue(:,:,:,1))     
      write(LDT_logunit,*) "Done reading HT_ROOF values."
      call read_netcdf_4d(n,trim(CLM45_struc(n)%surfacefile),"THICK_ROOF",3, &
          "numurbl","lsmlat","lsmlon","none", &
          LDT_rc%lnc(n),LDT_rc%lnr(n),CLM45_struc(n)%numurbl,1, &
          CLM45_struc(n)%THICK_ROOF%dvalue(:,:,:,1))     
      write(LDT_logunit,*) "Done reading THICK_ROOF values."
      call read_netcdf_4d(n,trim(CLM45_struc(n)%surfacefile),"THICK_WALL",3, &
          "numurbl","lsmlat","lsmlon","none", &
          LDT_rc%lnc(n),LDT_rc%lnr(n),CLM45_struc(n)%numurbl,1, &
          CLM45_struc(n)%THICK_WALL%dvalue(:,:,:,1))     
      write(LDT_logunit,*) "Done reading THICK_WALL values."
      call read_netcdf_4d(n,trim(CLM45_struc(n)%surfacefile),"T_BUILDING_MAX",3, &
          "numurbl","lsmlat","lsmlon","none", &
          LDT_rc%lnc(n),LDT_rc%lnr(n),CLM45_struc(n)%numurbl,1, &
          CLM45_struc(n)%T_BUILDING_MAX%dvalue(:,:,:,1))     
      write(LDT_logunit,*) "Done reading T_BUILDING_MAX values."
      call read_netcdf_4d(n,trim(CLM45_struc(n)%surfacefile),"T_BUILDING_MIN",3, &
          "numurbl","lsmlat","lsmlon","none", &
          LDT_rc%lnc(n),LDT_rc%lnr(n),CLM45_struc(n)%numurbl,1, &
          CLM45_struc(n)%T_BUILDING_MIN%dvalue(:,:,:,1))     
      write(LDT_logunit,*) "Done reading T_BUILDING_MIN values."
      call read_netcdf_4d(n,trim(CLM45_struc(n)%surfacefile),"WIND_HGT_CANYON",3, &
          "numurbl","lsmlat","lsmlon","none", &
          LDT_rc%lnc(n),LDT_rc%lnr(n),CLM45_struc(n)%numurbl,1, &
          CLM45_struc(n)%WIND_HGT_CANYON%dvalue(:,:,:,1))     
      write(LDT_logunit,*) "Done reading WIND_HGT_CANYON values."
      call read_netcdf_4d(n,trim(CLM45_struc(n)%surfacefile),"WTLUNIT_ROOF",3, &
          "numurbl","lsmlat","lsmlon","none", &
          LDT_rc%lnc(n),LDT_rc%lnr(n),CLM45_struc(n)%numurbl,1, &
          CLM45_struc(n)%WTLUNIT_ROOF%dvalue(:,:,:,1))     
      write(LDT_logunit,*) "Done reading WTLUNIT_ROOF values."
      call read_netcdf_4d(n,trim(CLM45_struc(n)%surfacefile),"WTROAD_PERV",3, &
          "numurbl","lsmlat","lsmlon","none", &
          LDT_rc%lnc(n),LDT_rc%lnr(n),CLM45_struc(n)%numurbl,1, &
          CLM45_struc(n)%WTROAD_PERV%dvalue(:,:,:,1))     
      write(LDT_logunit,*) "Done reading WTROAD_PERV values."
      call read_netcdf_4d(n,trim(CLM45_struc(n)%surfacefile),"ALB_IMPROAD_DIR",4, &
          "numrad","numurbl","lsmlat","lsmlon", &
          LDT_rc%lnc(n),LDT_rc%lnr(n),CLM45_struc(n)%numurbl,CLM45_struc(n)%numrad, &
          CLM45_struc(n)%ALB_IMPROAD_DIR%dvalue(:,:,:,:))     
      write(LDT_logunit,*) "Done reading ALB_IMPROAD_DIR values."
      call read_netcdf_4d(n,trim(CLM45_struc(n)%surfacefile),"ALB_IMPROAD_DIF",4, &
          "numrad","numurbl","lsmlat","lsmlon", &
          LDT_rc%lnc(n),LDT_rc%lnr(n),CLM45_struc(n)%numurbl,CLM45_struc(n)%numrad, &
          CLM45_struc(n)%ALB_IMPROAD_DIF%dvalue(:,:,:,:))     
      write(LDT_logunit,*) "Done reading ALB_IMPROAD_DIF values."
      call read_netcdf_4d(n,trim(CLM45_struc(n)%surfacefile),"ALB_PERROAD_DIR",4, &
          "numrad","numurbl","lsmlat","lsmlon", &
          LDT_rc%lnc(n),LDT_rc%lnr(n),CLM45_struc(n)%numurbl,CLM45_struc(n)%numrad, &
          CLM45_struc(n)%ALB_PERROAD_DIR%dvalue(:,:,:,:))     
      write(LDT_logunit,*) "Done reading ALB_PERROAD_DIR values."
      call read_netcdf_4d(n,trim(CLM45_struc(n)%surfacefile),"ALB_PERROAD_DIF",4, &
          "numrad","numurbl","lsmlat","lsmlon", &
          LDT_rc%lnc(n),LDT_rc%lnr(n),CLM45_struc(n)%numurbl,CLM45_struc(n)%numrad, &
          CLM45_struc(n)%ALB_PERROAD_DIF%dvalue(:,:,:,:))     
      write(LDT_logunit,*) "Done reading ALB_PERROAD_DIF values."
      call read_netcdf_4d(n,trim(CLM45_struc(n)%surfacefile),"ALB_ROOF_DIR",4, &
          "numrad","numurbl","lsmlat","lsmlon", &
          LDT_rc%lnc(n),LDT_rc%lnr(n),CLM45_struc(n)%numurbl,CLM45_struc(n)%numrad, &
          CLM45_struc(n)%ALB_ROOF_DIR%dvalue(:,:,:,:))     
      write(LDT_logunit,*) "Done reading ALB_ROOF_DIR values."
      call read_netcdf_4d(n,trim(CLM45_struc(n)%surfacefile),"ALB_ROOF_DIF",4, &
          "numrad","numurbl","lsmlat","lsmlon", &
          LDT_rc%lnc(n),LDT_rc%lnr(n),CLM45_struc(n)%numurbl,CLM45_struc(n)%numrad, &
          CLM45_struc(n)%ALB_ROOF_DIF%dvalue(:,:,:,:))     
      write(LDT_logunit,*) "Done reading ALB_ROOF_DIF values."
      call read_netcdf_4d(n,trim(CLM45_struc(n)%surfacefile),"ALB_WALL_DIR",4, &
          "numrad","numurbl","lsmlat","lsmlon", &
          LDT_rc%lnc(n),LDT_rc%lnr(n),CLM45_struc(n)%numurbl,CLM45_struc(n)%numrad, &
          CLM45_struc(n)%ALB_WALL_DIR%dvalue(:,:,:,:))     
      write(LDT_logunit,*) "Done reading ALB_WALL_DIR values."
      call read_netcdf_4d(n,trim(CLM45_struc(n)%surfacefile),"ALB_WALL_DIF",4, &
          "numrad","numurbl","lsmlat","lsmlon", &
          LDT_rc%lnc(n),LDT_rc%lnr(n),CLM45_struc(n)%numurbl,CLM45_struc(n)%numrad, &
          CLM45_struc(n)%ALB_WALL_DIF%dvalue(:,:,:,:))     
      write(LDT_logunit,*) "Done reading ALB_WALL_DIF values."
      call read_netcdf_4d(n,trim(CLM45_struc(n)%surfacefile),"TK_ROOF",4, &
          "nlevurb","numurbl","lsmlat","lsmlon", &
          LDT_rc%lnc(n),LDT_rc%lnr(n),CLM45_struc(n)%numurbl,CLM45_struc(n)%nlevurb, &
          CLM45_struc(n)%TK_ROOF%dvalue(:,:,:,:))     
      write(LDT_logunit,*) "Done reading TK_ROOF values."
      call read_netcdf_4d(n,trim(CLM45_struc(n)%surfacefile),"TK_WALL",4, &
          "nlevurb","numurbl","lsmlat","lsmlon", &
          LDT_rc%lnc(n),LDT_rc%lnr(n),CLM45_struc(n)%numurbl,CLM45_struc(n)%nlevurb, &
          CLM45_struc(n)%TK_WALL%dvalue(:,:,:,:))     
      write(LDT_logunit,*) "Done reading TK_WALL values."
      call read_netcdf_4d(n,trim(CLM45_struc(n)%surfacefile),"TK_IMPROAD",4, &
          "nlevurb","numurbl","lsmlat","lsmlon", &
          LDT_rc%lnc(n),LDT_rc%lnr(n),CLM45_struc(n)%numurbl,CLM45_struc(n)%nlevurb, &
          CLM45_struc(n)%TK_IMPROAD%dvalue(:,:,:,:))     
      write(LDT_logunit,*) "Done reading TK_IMPROAD values."
      call read_netcdf_4d(n,trim(CLM45_struc(n)%surfacefile),"CV_ROOF",4, &
          "nlevurb","numurbl","lsmlat","lsmlon", &
          LDT_rc%lnc(n),LDT_rc%lnr(n),CLM45_struc(n)%numurbl,CLM45_struc(n)%nlevurb, &
          CLM45_struc(n)%CV_ROOF%dvalue(:,:,:,:))     
      write(LDT_logunit,*) "Done reading CV_ROOF values."
      call read_netcdf_4d(n,trim(CLM45_struc(n)%surfacefile),"CV_WALL",4, &
          "nlevurb","numurbl","lsmlat","lsmlon", &
          LDT_rc%lnc(n),LDT_rc%lnr(n),CLM45_struc(n)%numurbl,CLM45_struc(n)%nlevurb, &
          CLM45_struc(n)%CV_WALL%dvalue(:,:,:,:))     
      write(LDT_logunit,*) "Done reading CV_WALL values."
      call read_netcdf_4d(n,trim(CLM45_struc(n)%surfacefile),"CV_IMPROAD",4, &
          "nlevurb","numurbl","lsmlat","lsmlon", &
          LDT_rc%lnc(n),LDT_rc%lnr(n),CLM45_struc(n)%numurbl,CLM45_struc(n)%nlevurb, &
          CLM45_struc(n)%CV_IMPROAD%dvalue(:,:,:,:))     
      write(LDT_logunit,*) "Done reading CV_IMPROAD values."
      call read_netcdf_4d(n,trim(CLM45_struc(n)%surfacefile),"NLEV_IMPROAD",3, &
          "numurbl","lsmlat","lsmlon","none", &
          LDT_rc%lnc(n),LDT_rc%lnr(n),CLM45_struc(n)%numurbl,1, &
          CLM45_struc(n)%NLEV_IMPROAD%dvalue(:,:,:,1))     
      write(LDT_logunit,*) "Done reading NLEV_IMPROAD values."
      call read_netcdf_4d(n,trim(CLM45_struc(n)%surfacefile),"peatf",2, &
          "lsmlat","lsmlon","none","none",LDT_rc%lnc(n),LDT_rc%lnr(n),1,1,&
          CLM45_struc(n)%peatf%dvalue(:,:,1,1))     
      write(LDT_logunit,*) "Done reading peatf values."
      call read_netcdf_4d(n,trim(CLM45_struc(n)%surfacefile),"abm",2, &
          "lsmlat","lsmlon","none","none",LDT_rc%lnc(n),LDT_rc%lnr(n),1,1,&
          CLM45_struc(n)%abm%dvalue(:,:,1,1))     
      write(LDT_logunit,*) "Done reading abm values."
      call read_netcdf_4d(n,trim(CLM45_struc(n)%surfacefile),"gdp",2, &
          "lsmlat","lsmlon","none","none",LDT_rc%lnc(n),LDT_rc%lnr(n),1,1,&
          CLM45_struc(n)%gdp%dvalue(:,:,1,1))     
      write(LDT_logunit,*) "Done reading gdp values."
      call read_netcdf_4d(n,trim(CLM45_struc(n)%surfacefile),"SLOPE",2, &
          "lsmlat","lsmlon","none","none",LDT_rc%lnc(n),LDT_rc%lnr(n),1,1,&
          CLM45_struc(n)%SLOPE%dvalue(:,:,1,1))     
      write(LDT_logunit,*) "Done reading SLOPE values."
      call read_netcdf_4d(n,trim(CLM45_struc(n)%surfacefile),"STD_ELEV",2, &
          "lsmlat","lsmlon","none","none",LDT_rc%lnc(n),LDT_rc%lnr(n),1,1,&
          CLM45_struc(n)%STD_ELEV%dvalue(:,:,1,1))     
      write(LDT_logunit,*) "Done reading STD_ELEV values."
      call read_netcdf_4d(n,trim(CLM45_struc(n)%surfacefile),"binfl",2, &
          "lsmlat","lsmlon","none","none",LDT_rc%lnc(n),LDT_rc%lnr(n),1,1,&
          CLM45_struc(n)%binfl%dvalue(:,:,1,1))     
      write(LDT_logunit,*) "Done reading binfl values."
      call read_netcdf_4d(n,trim(CLM45_struc(n)%surfacefile),"Ws",2, &
          "lsmlat","lsmlon","none","none",LDT_rc%lnc(n),LDT_rc%lnr(n),1,1,&
          CLM45_struc(n)%Ws%dvalue(:,:,1,1))     
      write(LDT_logunit,*) "Done reading Ws values."
      call read_netcdf_4d(n,trim(CLM45_struc(n)%surfacefile),"Dsmax",2, &
          "lsmlat","lsmlon","none","none",LDT_rc%lnc(n),LDT_rc%lnr(n),1,1,&
          CLM45_struc(n)%Dsmax%dvalue(:,:,1,1))     
      write(LDT_logunit,*) "Done reading Dsmax values."
      call read_netcdf_4d(n,trim(CLM45_struc(n)%surfacefile),"Ds",2, &
          "lsmlat","lsmlon","none","none",LDT_rc%lnc(n),LDT_rc%lnr(n),1,1,&
          CLM45_struc(n)%Ds%dvalue(:,:,1,1))     
      write(LDT_logunit,*) "Done reading Ds values."
      call read_netcdf_4d(n,trim(CLM45_struc(n)%surfacefile),"LAKEDEPTH",2, &
          "lsmlat","lsmlon","none","none",LDT_rc%lnc(n),LDT_rc%lnr(n),1,1,&
          CLM45_struc(n)%LAKEDEPTH%dvalue(:,:,1,1))     
      write(LDT_logunit,*) "Done reading LAKEDEPTH values."
      call read_netcdf_4d(n,trim(CLM45_struc(n)%surfacefile),"F0",2, &
          "lsmlat","lsmlon","none","none",LDT_rc%lnc(n),LDT_rc%lnr(n),1,1,&
          CLM45_struc(n)%F0%dvalue(:,:,1,1))     
      write(LDT_logunit,*) "Done reading F0 values."
      call read_netcdf_4d(n,trim(CLM45_struc(n)%surfacefile),"P3",2, &
          "lsmlat","lsmlon","none","none",LDT_rc%lnc(n),LDT_rc%lnr(n),1,1,&
          CLM45_struc(n)%P3%dvalue(:,:,1,1))     
      write(LDT_logunit,*) "Done reading P3 values."
      call read_netcdf_4d(n,trim(CLM45_struc(n)%surfacefile),"ZWT0",2, &
          "lsmlat","lsmlon","none","none",LDT_rc%lnc(n),LDT_rc%lnr(n),1,1,&
          CLM45_struc(n)%ZWT0%dvalue(:,:,1,1))     
      write(LDT_logunit,*) "Done reading ZWT0 values."
      call read_netcdf_4d(n,trim(CLM45_struc(n)%surfacefile),"PCT_WETLAND",2, &
          "lsmlat","lsmlon","none","none",LDT_rc%lnc(n),LDT_rc%lnr(n),1,1,&
          CLM45_struc(n)%PCT_WETLAND%dvalue(:,:,1,1))     
      write(LDT_logunit,*) "Done reading PCT_WETLAND values."
      call read_netcdf_4d(n,trim(CLM45_struc(n)%surfacefile),"PCT_LAKE",2, &
          "lsmlat","lsmlon","none","none",LDT_rc%lnc(n),LDT_rc%lnr(n),1,1,&
          CLM45_struc(n)%PCT_LAKE%dvalue(:,:,1,1))     
      write(LDT_logunit,*) "Done reading PCT_LAKE values."
      call read_netcdf_4d(n,trim(CLM45_struc(n)%surfacefile),"PCT_GLACIER",2, &
          "lsmlat","lsmlon","none","none",LDT_rc%lnc(n),LDT_rc%lnr(n),1,1,&
          CLM45_struc(n)%PCT_GLACIER%dvalue(:,:,1,1))     
      write(LDT_logunit,*) "Done reading PCT_GLACIER values."
      call read_netcdf_4d(n,trim(CLM45_struc(n)%surfacefile),"GLC_MEC",1, &
          "nglcecp1","none","none","none",CLM45_struc(n)%nglcecp1,1,1,1,&
          CLM45_struc(n)%GLC_MEC%dvalue(:,1,1,1))     
      write(LDT_logunit,*) "Done reading GLC_MEC values."
      call read_netcdf_4d(n,trim(CLM45_struc(n)%surfacefile),"PCT_GLC_MEC",3, &
          "nglcec","lsmlat","lsmlon","none", &
          LDT_rc%lnc(n),LDT_rc%lnr(n),CLM45_struc(n)%nglcec,1, &
          CLM45_struc(n)%PCT_GLC_MEC%dvalue(:,:,:,1))     
      write(LDT_logunit,*) "Done reading PCT_GLC_MEC values."
      call read_netcdf_4d(n,trim(CLM45_struc(n)%surfacefile),"PCT_GLC_MEC_GIC",3, &
          "nglcec","lsmlat","lsmlon","none", &
          LDT_rc%lnc(n),LDT_rc%lnr(n),CLM45_struc(n)%nglcec,1, &
          CLM45_struc(n)%PCT_GLC_MEC_GIC%dvalue(:,:,:,1))     
      write(LDT_logunit,*) "Done reading PCT_GLC_MEC_GIC values."
      call read_netcdf_4d(n,trim(CLM45_struc(n)%surfacefile),"PCT_GLC_MEC_ICESHEET",3, &
          "nglcec","lsmlat","lsmlon","none", &
          LDT_rc%lnc(n),LDT_rc%lnr(n),CLM45_struc(n)%nglcec,1, &
          CLM45_struc(n)%PCT_GLC_MEC_ICESHEET%dvalue(:,:,:,1))     
      write(LDT_logunit,*) "Done reading PCT_GLC_MEC_ICESHEET values."
      call read_netcdf_4d(n,trim(CLM45_struc(n)%surfacefile),"PCT_GLC_GIC",2, &
          "lsmlat","lsmlon","none","none",LDT_rc%lnc(n),LDT_rc%lnr(n),1,1,&
          CLM45_struc(n)%PCT_GLC_GIC%dvalue(:,:,1,1))     
      write(LDT_logunit,*) "Done reading PCT_GLC_GIC values."
      call read_netcdf_4d(n,trim(CLM45_struc(n)%surfacefile),"PCT_GLC_ICESHEET",2, &
          "lsmlat","lsmlon","none","none",LDT_rc%lnc(n),LDT_rc%lnr(n),1,1,&
          CLM45_struc(n)%PCT_GLC_ICESHEET%dvalue(:,:,1,1))     
      write(LDT_logunit,*) "Done reading PCT_GLC_ICESHEET values."
      call read_netcdf_4d(n,trim(CLM45_struc(n)%surfacefile),"TOPO_GLC_MEC",3, &
          "nglcec","lsmlat","lsmlon","none", &
          LDT_rc%lnc(n),LDT_rc%lnr(n),CLM45_struc(n)%nglcec,1, &
          CLM45_struc(n)%TOPO_GLC_MEC%dvalue(:,:,:,1))     
      write(LDT_logunit,*) "Done reading TOPO_GLC_MEC values."
      call read_netcdf_4d(n,trim(CLM45_struc(n)%surfacefile),"TOPO",2, &
          "lsmlat","lsmlon","none","none",LDT_rc%lnc(n),LDT_rc%lnr(n),1,1,&
          CLM45_struc(n)%TOPO%dvalue(:,:,1,1))     
      write(LDT_logunit,*) "Done reading TOPO values."
      call read_netcdf_4d(n,trim(CLM45_struc(n)%surfacefile),"PCT_URBAN",3, &
          "numurbl","lsmlat","lsmlon","none", &
          LDT_rc%lnc(n),LDT_rc%lnr(n),CLM45_struc(n)%numurbl,1, &
          CLM45_struc(n)%PCT_URBAN%dvalue(:,:,:,1))     
      write(LDT_logunit,*) "Done reading PCT_URBAN values."
      call read_netcdf_4d(n,trim(CLM45_struc(n)%surfacefile),"URBAN_REGION_ID",2, &
          "lsmlat","lsmlon","none","none",LDT_rc%lnc(n),LDT_rc%lnr(n),1,1,&
          CLM45_struc(n)%URBAN_REGION_ID%dvalue(:,:,1,1))     
      write(LDT_logunit,*) "Done reading URBAN_REGION_ID values."

 end subroutine clm45_read_params
 
 subroutine getclm45lcst(n)
! assign grid% for land cover fraction and set surfacetype
! replaces Yudong's single_col.F90 routine
   integer   :: n 
   integer :: surfacetype(LDT_rc%lnc(n),LDT_rc%lnr(n),LDT_LSMparam_struc(n)%landcover%num_bins)
   real*8  :: landcover(LDT_rc%lnc(n),LDT_rc%lnr(n),LDT_LSMparam_struc(n)%landcover%num_bins)   ! grid fractions of each pft, lake, wetland, & glacier
   integer :: ic, ir, iz, t

   surfacetype = 0
   landcover = 0.0

   do ir = 1, LDT_rc%lnr(n)
    do ic = 1, LDT_rc%lnc(n)
      !if( LDT_rc%global_mask(ic,ir) >= 1 ) then
      !==> use CLM-4.5 mask already cropped for subdomain
      if( CLM45_struc(n)%mask%dvalue(ic,ir,1,1) >= 1 ) then
        !vegetated set all surfacetype to be 1
        if ( sum(CLM45_struc(n)%PCT_PFT%dvalue(ic,ir,1:17,1)) .gt. 0.0 ) then
          surfacetype(ic, ir, 1:17) = 1
        endif
        do iz = 1, 17  ! pft 1-17
          if ( CLM45_struc(n)%PCT_PFT%dvalue(ic,ir,iz,1) .gt. 0.0 ) then
           landcover(ic, ir, iz) = CLM45_struc(n)%PCT_PFT%dvalue(ic,ir,iz,1) *0.01  !in fraction
          endif
        end do   ! iz
        !urban_tbd iz=18:22
        if ( CLM45_struc(n)%PCT_URBAN%dvalue(ic,ir,1,1) .gt. 0.0 ) then
           landcover(ic, ir, 18) = 0.01*CLM45_struc(n)%PCT_URBAN%dvalue(ic,ir,1,1)*CLM45_struc(n)%WTLUNIT_ROOF%dvalue(ic,ir,1,1)
           landcover(ic, ir, 19) = 0.01*CLM45_struc(n)%PCT_URBAN%dvalue(ic,ir,1,1)*(1.0-CLM45_struc(n)%WTLUNIT_ROOF%dvalue(ic,ir,1,1))/3.0 
           landcover(ic, ir, 20) = 0.01*CLM45_struc(n)%PCT_URBAN%dvalue(ic,ir,1,1)*(1.0-CLM45_struc(n)%WTLUNIT_ROOF%dvalue(ic,ir,1,1))/3.0 
           landcover(ic, ir, 21) = 0.01*CLM45_struc(n)%PCT_URBAN%dvalue(ic,ir,1,1)*(1.0-CLM45_struc(n)%WTLUNIT_ROOF%dvalue(ic,ir,1,1))/3.0*(1.0-CLM45_struc(n)%WTROAD_PERV%dvalue(ic,ir,1,1)) 
           landcover(ic, ir, 22) = 0.01*CLM45_struc(n)%PCT_URBAN%dvalue(ic,ir,1,1)*(1.0-CLM45_struc(n)%WTLUNIT_ROOF%dvalue(ic,ir,1,1))/3.0*CLM45_struc(n)%WTROAD_PERV%dvalue(ic,ir,1,1)
           if (landcover(ic, ir, 18).gt.0.0) surfacetype(ic, ir, 18) = 1
           if (landcover(ic, ir, 19).gt.0.0) surfacetype(ic, ir, 19) = 1
           if (landcover(ic, ir, 20).gt.0.0) surfacetype(ic, ir, 20) = 1
           if (landcover(ic, ir, 21).gt.0.0) surfacetype(ic, ir, 21) = 1
           if (landcover(ic, ir, 22).gt.0.0) surfacetype(ic, ir, 22) = 1
        endif
        !urban_hd iz=23:27
        if ( CLM45_struc(n)%PCT_URBAN%dvalue(ic,ir,2,1) .gt. 0.0 ) then
           landcover(ic, ir, 23) = 0.01*CLM45_struc(n)%PCT_URBAN%dvalue(ic,ir,2,1)*CLM45_struc(n)%WTLUNIT_ROOF%dvalue(ic,ir,2,1)
           landcover(ic, ir, 24) = 0.01*CLM45_struc(n)%PCT_URBAN%dvalue(ic,ir,2,1)*(1.0-CLM45_struc(n)%WTLUNIT_ROOF%dvalue(ic,ir,2,1))/3.0 
           landcover(ic, ir, 25) = 0.01*CLM45_struc(n)%PCT_URBAN%dvalue(ic,ir,2,1)*(1.0-CLM45_struc(n)%WTLUNIT_ROOF%dvalue(ic,ir,2,1))/3.0 
           landcover(ic, ir, 26) = 0.01*CLM45_struc(n)%PCT_URBAN%dvalue(ic,ir,2,1)*(1.0-CLM45_struc(n)%WTLUNIT_ROOF%dvalue(ic,ir,2,1))/3.0*(1.0-CLM45_struc(n)%WTROAD_PERV%dvalue(ic,ir,2,1)) 
           landcover(ic, ir, 27) = 0.01*CLM45_struc(n)%PCT_URBAN%dvalue(ic,ir,2,1)*(1.0-CLM45_struc(n)%WTLUNIT_ROOF%dvalue(ic,ir,2,1))/3.0*CLM45_struc(n)%WTROAD_PERV%dvalue(ic,ir,2,1)
           if (landcover(ic, ir, 23).gt.0.0) surfacetype(ic, ir, 23) = 1
           if (landcover(ic, ir, 24).gt.0.0) surfacetype(ic, ir, 24) = 1
           if (landcover(ic, ir, 25).gt.0.0) surfacetype(ic, ir, 25) = 1
           if (landcover(ic, ir, 26).gt.0.0) surfacetype(ic, ir, 26) = 1
           if (landcover(ic, ir, 27).gt.0.0) surfacetype(ic, ir, 27) = 1
        endif
        !urban_md iz=28:32
        if ( CLM45_struc(n)%PCT_URBAN%dvalue(ic,ir,3,1) .gt. 0.0 ) then
           landcover(ic, ir, 28) = 0.01*CLM45_struc(n)%PCT_URBAN%dvalue(ic,ir,3,1)*CLM45_struc(n)%WTLUNIT_ROOF%dvalue(ic,ir,3,1)
           landcover(ic, ir, 29) = 0.01*CLM45_struc(n)%PCT_URBAN%dvalue(ic,ir,3,1)*(1.0-CLM45_struc(n)%WTLUNIT_ROOF%dvalue(ic,ir,3,1))/3.0 
           landcover(ic, ir, 30) = 0.01*CLM45_struc(n)%PCT_URBAN%dvalue(ic,ir,3,1)*(1.0-CLM45_struc(n)%WTLUNIT_ROOF%dvalue(ic,ir,3,1))/3.0 
           landcover(ic, ir, 31) = 0.01*CLM45_struc(n)%PCT_URBAN%dvalue(ic,ir,3,1)*(1.0-CLM45_struc(n)%WTLUNIT_ROOF%dvalue(ic,ir,3,1))/3.0*(1.0-CLM45_struc(n)%WTROAD_PERV%dvalue(ic,ir,3,1)) 
           landcover(ic, ir, 32) = 0.01*CLM45_struc(n)%PCT_URBAN%dvalue(ic,ir,3,1)*(1.0-CLM45_struc(n)%WTLUNIT_ROOF%dvalue(ic,ir,3,1))/3.0*CLM45_struc(n)%WTROAD_PERV%dvalue(ic,ir,3,1)
           if (landcover(ic, ir, 28).gt.0.0) surfacetype(ic, ir, 28) = 1
           if (landcover(ic, ir, 29).gt.0.0) surfacetype(ic, ir, 29) = 1
           if (landcover(ic, ir, 30).gt.0.0) surfacetype(ic, ir, 30) = 1
           if (landcover(ic, ir, 31).gt.0.0) surfacetype(ic, ir, 31) = 1
           if (landcover(ic, ir, 32).gt.0.0) surfacetype(ic, ir, 32) = 1
        endif
        ! lake iz=33
        if ( CLM45_struc(n)%PCT_LAKE%dvalue(ic,ir,1,1) .gt. 0.0 ) then
           surfacetype(ic, ir, 33) = 1
           landcover(ic, ir, 33) = CLM45_struc(n)%PCT_LAKE%dvalue(ic,ir,1,1)*0.01
        endif
        ! wetland iz=34
        if ( CLM45_struc(n)%PCT_WETLAND%dvalue(ic,ir,1,1) .gt. 0.0 ) then
           surfacetype(ic, ir, 34) = 1
           landcover(ic, ir, 34) = CLM45_struc(n)%PCT_WETLAND%dvalue(ic,ir,1,1)*0.01
        endif
        ! glacier iz=35
        if ( CLM45_struc(n)%PCT_GLACIER%dvalue(ic,ir,1,1) .gt. 0.0 ) then
           surfacetype(ic, ir, 35) = 1
           landcover(ic, ir, 35) = CLM45_struc(n)%PCT_GLACIER%dvalue(ic,ir,1,1)*0.01
        endif
        ! water iz=36
        surfacetype(ic, ir, 36) = 0
        landcover(ic, ir, 36) = 0.0
      endif   ! land
    end do  ! ic
   end do  ! ir

   do ir = 1, LDT_rc%lnr(n)
    do ic = 1, LDT_rc%lnc(n)
     do t = 1, LDT_LSMparam_struc(n)%landcover%num_bins   ! including 36=water
       LDT_LSMparam_struc(n)%landcover%value(ic,ir,t) = real(landcover(ic,ir,t))
       LDT_LSMparam_struc(n)%sfctype%value(ic,ir,t) = surfacetype(ic,ir,t)*1.0
     end do  ! t
    end do  ! ic
   end do  ! ir
  
   write(LDT_logunit,*)" ** MSG:  done getclm45lcst"
    
 end subroutine getclm45lcst

 subroutine CLM45parms_writeHeader(n,ftn,dimID,monthID)
   
   integer   :: n 
   integer   :: ftn
   integer   :: dimID(3)
   integer   :: monthID
   
   integer   :: ndimID(3)  ! 3D, vlevel>1
   integer   :: fdimID(4)  ! 4D, vlevel>1 & zlevel>1
   integer   :: tempID, udimID
   
   if( CLM45_struc(n)%clm45%selectOpt == 1 ) then
    fdimID(1) = dimID(1)
    fdimID(2) = dimID(2)
    fdimID(4) = monthID


    call LDT_writeNETCDFdataHeader(n,ftn,dimID,&
         CLM45_struc(n)%xc,"double")
    call LDT_writeNETCDFdataHeader(n,ftn,dimID,&
         CLM45_struc(n)%yc,"double")
    call LDT_verify(nf90_def_dim(ftn,'ncorners',&
          CLM45_struc(n)%ncorners,ndimID(3)))
    ndimID(1) = dimID(1)
    ndimID(2) = dimID(2)
    call LDT_writeNETCDFdataHeader(n,ftn,ndimID,&
         CLM45_struc(n)%xv,"double")
    call LDT_writeNETCDFdataHeader(n,ftn,ndimID,&
         CLM45_struc(n)%yv,"double")
    call LDT_writeNETCDFdataHeader(n,ftn,dimID,&
         CLM45_struc(n)%mask,"int")
    call LDT_writeNETCDFdataHeader(n,ftn,dimID,&
         CLM45_struc(n)%area,"double")
    call LDT_writeNETCDFdataHeader(n,ftn,dimID,&
         CLM45_struc(n)%frac,"double")
    call LDT_verify(nf90_def_var(ftn, &
         trim(CLM45_struc(n)%mxsoil_color%short_name), &
         nf90_int,varID=CLM45_struc(n)%mxsoil_color%vid))
    call LDT_writeNETCDFdataHeader(n,ftn,dimID,&
         CLM45_struc(n)%SOIL_COLOR,"int")
    call LDT_verify(nf90_def_dim(ftn,'nlevsoi',&
          CLM45_struc(n)%nlevsoi,ndimID(3)))
    ndimID(1) = dimID(1)
    ndimID(2) = dimID(2)
    call LDT_writeNETCDFdataHeader(n,ftn,ndimID,&
         CLM45_struc(n)%PCT_SAND,"double")
    call LDT_writeNETCDFdataHeader(n,ftn,ndimID,&
         CLM45_struc(n)%PCT_CLAY,"double")
    call LDT_writeNETCDFdataHeader(n,ftn,ndimID,&
         CLM45_struc(n)%ORGANIC,"double")
    call LDT_writeNETCDFdataHeader(n,ftn,dimID,&
         CLM45_struc(n)%FMAX,"double")
    call LDT_writeNETCDFdataHeader(n,ftn,dimID,&
         CLM45_struc(n)%LANDFRAC_PFT,"double")
    call LDT_writeNETCDFdataHeader(n,ftn,dimID,&
         CLM45_struc(n)%PFTDATA_MASK,"int")
    call LDT_verify(nf90_def_dim(ftn,'lsmpft',&
          CLM45_struc(n)%lsmpft,ndimID(3)))
    call LDT_writeNETCDFdataHeader(n,ftn,ndimID,&
         CLM45_struc(n)%PCT_PFT,"double")
    fdimID(3) = ndimID(3)   ! lsmpft
    call LDT_writeNETCDFdataHeader(n,ftn,fdimID(1:3),&
         CLM45_struc(n)%MONTHLY_LAI,"double",fdimID)
    call LDT_writeNETCDFdataHeader(n,ftn,fdimID(1:3),&
         CLM45_struc(n)%MONTHLY_SAI,"double",fdimID)
    call LDT_writeNETCDFdataHeader(n,ftn,fdimID(1:3),&
         CLM45_struc(n)%MONTHLY_HEIGHT_TOP,"double",fdimID)
    call LDT_writeNETCDFdataHeader(n,ftn,fdimID(1:3),&
         CLM45_struc(n)%MONTHLY_HEIGHT_BOT,"double",fdimID)
    ! skip int time(time), taken care of at higher level
    call LDT_writeNETCDFdataHeader(n,ftn,dimID,&
         CLM45_struc(n)%SRFAREA,"double")
    call LDT_writeNETCDFdataHeader(n,ftn,dimID,&
         CLM45_struc(n)%LONGXY,"double")
    call LDT_writeNETCDFdataHeader(n,ftn,dimID,&
         CLM45_struc(n)%LATIXY,"double")
    call LDT_writeNETCDFdataHeader(n,ftn,dimID,&
         CLM45_struc(n)%EF1_BTR,"double")
    call LDT_writeNETCDFdataHeader(n,ftn,dimID,&
         CLM45_struc(n)%EF1_FET,"double")
    call LDT_writeNETCDFdataHeader(n,ftn,dimID,&
         CLM45_struc(n)%EF1_FDT,"double")
    call LDT_writeNETCDFdataHeader(n,ftn,dimID,&
         CLM45_struc(n)%EF1_SHR,"double")
    call LDT_writeNETCDFdataHeader(n,ftn,dimID,&
         CLM45_struc(n)%EF1_GRS,"double")
    call LDT_writeNETCDFdataHeader(n,ftn,dimID,&
         CLM45_struc(n)%EF1_CRP,"double")
    call LDT_verify(nf90_def_dim(ftn,'numurbl',&
          CLM45_struc(n)%numurbl,ndimID(3)))
    call LDT_writeNETCDFdataHeader(n,ftn,ndimID,&
         CLM45_struc(n)%CANYON_HWR,"double")
    udimID = ndimID(3)  ! save numurbl
    call LDT_writeNETCDFdataHeader(n,ftn,ndimID,&
         CLM45_struc(n)%EM_IMPROAD,"double")
    call LDT_writeNETCDFdataHeader(n,ftn,ndimID,&
         CLM45_struc(n)%EM_PERROAD,"double")
    call LDT_writeNETCDFdataHeader(n,ftn,ndimID,&
         CLM45_struc(n)%EM_ROOF,"double")
    call LDT_writeNETCDFdataHeader(n,ftn,ndimID,&
         CLM45_struc(n)%EM_WALL,"double")
    call LDT_writeNETCDFdataHeader(n,ftn,ndimID,&
         CLM45_struc(n)%HT_ROOF,"double")
    call LDT_writeNETCDFdataHeader(n,ftn,ndimID,&
         CLM45_struc(n)%THICK_ROOF,"double")
    call LDT_writeNETCDFdataHeader(n,ftn,ndimID,&
         CLM45_struc(n)%THICK_WALL,"double")
    call LDT_writeNETCDFdataHeader(n,ftn,ndimID,&
         CLM45_struc(n)%T_BUILDING_MAX,"double")
    call LDT_writeNETCDFdataHeader(n,ftn,ndimID,&
         CLM45_struc(n)%T_BUILDING_MIN,"double")
    call LDT_writeNETCDFdataHeader(n,ftn,ndimID,&
         CLM45_struc(n)%WIND_HGT_CANYON,"double")
    call LDT_writeNETCDFdataHeader(n,ftn,ndimID,&
         CLM45_struc(n)%WTLUNIT_ROOF,"double")
    call LDT_writeNETCDFdataHeader(n,ftn,ndimID,&
         CLM45_struc(n)%WTROAD_PERV,"double")
    fdimID(3) = udimID   !numurbl
    call LDT_verify(nf90_def_dim(ftn,'numrad',&
          CLM45_struc(n)%numrad,fdimID(4)))
    call LDT_writeNETCDFdataHeader(n,ftn,fdimID(1:3),&
         CLM45_struc(n)%ALB_IMPROAD_DIR,"double",fdimID)
    call LDT_writeNETCDFdataHeader(n,ftn,fdimID(1:3),&
         CLM45_struc(n)%ALB_IMPROAD_DIF,"double",fdimID)
    call LDT_writeNETCDFdataHeader(n,ftn,fdimID(1:3),&
         CLM45_struc(n)%ALB_PERROAD_DIR,"double",fdimID)
    call LDT_writeNETCDFdataHeader(n,ftn,fdimID(1:3),&
         CLM45_struc(n)%ALB_PERROAD_DIF,"double",fdimID)
    call LDT_writeNETCDFdataHeader(n,ftn,fdimID(1:3),&
         CLM45_struc(n)%ALB_ROOF_DIR,"double",fdimID)
    call LDT_writeNETCDFdataHeader(n,ftn,fdimID(1:3),&
         CLM45_struc(n)%ALB_ROOF_DIF,"double",fdimID)
    call LDT_writeNETCDFdataHeader(n,ftn,fdimID(1:3),&
         CLM45_struc(n)%ALB_WALL_DIR,"double",fdimID)
    call LDT_writeNETCDFdataHeader(n,ftn,fdimID(1:3),&
         CLM45_struc(n)%ALB_WALL_DIF,"double",fdimID)
    call LDT_verify(nf90_def_dim(ftn,'nlevurb',&
          CLM45_struc(n)%nlevurb,fdimID(4)))
    call LDT_writeNETCDFdataHeader(n,ftn,fdimID(1:3),&
         CLM45_struc(n)%TK_ROOF,"double",fdimID)
    call LDT_writeNETCDFdataHeader(n,ftn,fdimID(1:3),&
         CLM45_struc(n)%TK_WALL,"double",fdimID)
    call LDT_writeNETCDFdataHeader(n,ftn,fdimID(1:3),&
         CLM45_struc(n)%TK_IMPROAD,"double",fdimID)
    call LDT_writeNETCDFdataHeader(n,ftn,fdimID(1:3),&
         CLM45_struc(n)%CV_ROOF,"double",fdimID)
    call LDT_writeNETCDFdataHeader(n,ftn,fdimID(1:3),&
         CLM45_struc(n)%CV_WALL,"double",fdimID)
    call LDT_writeNETCDFdataHeader(n,ftn,fdimID(1:3),&
         CLM45_struc(n)%CV_IMPROAD,"double",fdimID)
    call LDT_writeNETCDFdataHeader(n,ftn,ndimID,&
         CLM45_struc(n)%NLEV_IMPROAD,"int")
    call LDT_writeNETCDFdataHeader(n,ftn,dimID,&
         CLM45_struc(n)%peatf,"double")
    call LDT_writeNETCDFdataHeader(n,ftn,dimID,&
         CLM45_struc(n)%abm,"int")
    call LDT_writeNETCDFdataHeader(n,ftn,dimID,&
         CLM45_struc(n)%gdp,"double")
    call LDT_writeNETCDFdataHeader(n,ftn,dimID,&
         CLM45_struc(n)%SLOPE,"double")
    call LDT_writeNETCDFdataHeader(n,ftn,dimID,&
         CLM45_struc(n)%STD_ELEV,"double")
    call LDT_writeNETCDFdataHeader(n,ftn,dimID,&
         CLM45_struc(n)%binfl,"double")
    call LDT_writeNETCDFdataHeader(n,ftn,dimID,&
         CLM45_struc(n)%Ws,"double")
    call LDT_writeNETCDFdataHeader(n,ftn,dimID,&
         CLM45_struc(n)%Dsmax,"double")
    call LDT_writeNETCDFdataHeader(n,ftn,dimID,&
         CLM45_struc(n)%Ds,"double")
    call LDT_writeNETCDFdataHeader(n,ftn,dimID,&
         CLM45_struc(n)%LAKEDEPTH,"double")
    call LDT_writeNETCDFdataHeader(n,ftn,dimID,&
         CLM45_struc(n)%F0,"double")
    call LDT_writeNETCDFdataHeader(n,ftn,dimID,&
         CLM45_struc(n)%P3,"double")
    call LDT_writeNETCDFdataHeader(n,ftn,dimID,&
         CLM45_struc(n)%ZWT0,"double")
    call LDT_writeNETCDFdataHeader(n,ftn,dimID,&
         CLM45_struc(n)%PCT_WETLAND,"double")
    call LDT_writeNETCDFdataHeader(n,ftn,dimID,&
         CLM45_struc(n)%PCT_LAKE,"double")
    call LDT_writeNETCDFdataHeader(n,ftn,dimID,&
         CLM45_struc(n)%PCT_GLACIER,"double")
    call LDT_verify(nf90_def_dim(ftn,'nglcecp1',&
          CLM45_struc(n)%nglcecp1,tempID))
    call LDT_verify(nf90_def_var(ftn, &
         trim(CLM45_struc(n)%GLC_MEC%short_name), &
         nf90_double,dimids=tempID,varID=CLM45_struc(n)%GLC_MEC%vid))
    call LDT_verify(nf90_def_dim(ftn,'nglcec',&
          CLM45_struc(n)%nglcec,ndimID(3)))
    call LDT_writeNETCDFdataHeader(n,ftn,ndimID,&
         CLM45_struc(n)%PCT_GLC_MEC,"double")
    call LDT_writeNETCDFdataHeader(n,ftn,ndimID,&
         CLM45_struc(n)%PCT_GLC_MEC_GIC,"double")
    call LDT_writeNETCDFdataHeader(n,ftn,ndimID,&
         CLM45_struc(n)%PCT_GLC_MEC_ICESHEET,"double")
    call LDT_writeNETCDFdataHeader(n,ftn,dimID,&
         CLM45_struc(n)%PCT_GLC_GIC,"double")
    call LDT_writeNETCDFdataHeader(n,ftn,dimID,&
         CLM45_struc(n)%PCT_GLC_ICESHEET,"double")
    call LDT_writeNETCDFdataHeader(n,ftn,ndimID,&
         CLM45_struc(n)%TOPO_GLC_MEC,"double")
    call LDT_writeNETCDFdataHeader(n,ftn,dimID,&
         CLM45_struc(n)%TOPO,"double")
    ndimID(3) = udimID   ! numurbl
    call LDT_writeNETCDFdataHeader(n,ftn,ndimID,&
         CLM45_struc(n)%PCT_URBAN,"double")
    call LDT_writeNETCDFdataHeader(n,ftn,dimID,&
         CLM45_struc(n)%URBAN_REGION_ID,"int")
!    print*,'done writeNETCDFdataHeader'

   endif
   
 end subroutine CLM45parms_writeHeader
 
  subroutine CLM45parms_writeData(n,ftn)

    integer   :: n 
    integer   :: ftn

    if( CLM45_struc(n)%clm45%selectOpt == 1 ) then
      call LDT_writeNETCDFdata(n,ftn,CLM45_struc(n)%xc,"double")
!      print*,'wrote xc'
      call LDT_writeNETCDFdata(n,ftn,CLM45_struc(n)%yc,"double")
!      print*,'wrote yc'
      call LDT_writeNETCDFdata(n,ftn,CLM45_struc(n)%xv,"double")
!      print*,'wrote xv'
      call LDT_writeNETCDFdata(n,ftn,CLM45_struc(n)%yv,"double")
!      print*,'wrote yv'
      call LDT_writeNETCDFdata(n,ftn,CLM45_struc(n)%mask,"int")
!      print*,'wrote mask'
      call LDT_writeNETCDFdata(n,ftn,CLM45_struc(n)%area,"double")
!      print*,'wrote area'
      call LDT_writeNETCDFdata(n,ftn,CLM45_struc(n)%frac,"double")
!      print*,'wrote frac'
      ! special single variable: mxsoil_color
      call LDT_verify(nf90_put_var(ftn, CLM45_struc(n)%mxsoil_color%vid,&
           CLM45_struc(n)%mxsoil_color%dvalue), &
           'nf90_put_var failed in CLM45parms_writeData')
         call LDT_verify(nf90_put_att(ftn,CLM45_struc(n)%mxsoil_color%vid,&
              "standard_name",trim(CLM45_struc(n)%mxsoil_color%standard_name)))
         call LDT_verify(nf90_put_att(ftn,CLM45_struc(n)%mxsoil_color%vid,&
              "units",trim(CLM45_struc(n)%mxsoil_color%units)))
!      print*,'wrote mxsoil_color'
      call LDT_writeNETCDFdata(n,ftn,CLM45_struc(n)%SOIL_COLOR,"int")
!      print*,'wrote SOIL_COLOR'
      call LDT_writeNETCDFdata(n,ftn,CLM45_struc(n)%PCT_SAND,"double")
!      print*,'wrote PCT_SAND'
      call LDT_writeNETCDFdata(n,ftn,CLM45_struc(n)%PCT_CLAY,"double")
!      print*,'wrote PCT_CLAY'
      call LDT_writeNETCDFdata(n,ftn,CLM45_struc(n)%ORGANIC,"double")
!      print*,'wrote ORGANIC'
      call LDT_writeNETCDFdata(n,ftn,CLM45_struc(n)%FMAX,"double")
!      print*,'wrote FMAX'
      call LDT_writeNETCDFdata(n,ftn,CLM45_struc(n)%LANDFRAC_PFT,&
                               "double")
!      print*,'wrote LANDFRAC_PFT'
      call LDT_writeNETCDFdata(n,ftn,CLM45_struc(n)%PFTDATA_MASK,"int")
!      print*,'wrote PFTDATA_MASK'
      call LDT_writeNETCDFdata(n,ftn,CLM45_struc(n)%PCT_PFT,"double")
!      print*,'wrote PCT_PFT'
      call LDT_writeNETCDFdata(n,ftn,CLM45_struc(n)%MONTHLY_LAI,"double")
!      print*,'wrote MONTHLY_LAI'
      call LDT_writeNETCDFdata(n,ftn,CLM45_struc(n)%MONTHLY_SAI,"double")
!      print*,'wrote MONTHLY_SAI'
      call LDT_writeNETCDFdata(n,ftn,CLM45_struc(n)%MONTHLY_HEIGHT_TOP,"double")
!      print*,'wrote MONTHLY_HEIGHT_TOP'
      call LDT_writeNETCDFdata(n,ftn,CLM45_struc(n)%MONTHLY_HEIGHT_BOT,"double")
!      print*,'wrote MONTHLY_HEIGHT_BOT'
      ! skip time <- taken care of in historyMod
      call LDT_writeNETCDFdata(n,ftn,CLM45_struc(n)%SRFAREA,"double")
!      print*,'wrote SRFAREA'
      call LDT_writeNETCDFdata(n,ftn,CLM45_struc(n)%LONGXY,"double")
!      print*,'wrote LONGXY'
      call LDT_writeNETCDFdata(n,ftn,CLM45_struc(n)%LATIXY,"double")
!      print*,'wrote LATIXY'
      call LDT_writeNETCDFdata(n,ftn,CLM45_struc(n)%EF1_BTR,"double")
!      print*,'wrote EF1_BTR'
      call LDT_writeNETCDFdata(n,ftn,CLM45_struc(n)%EF1_FET,"double")
!      print*,'wrote EF1_FET'
      call LDT_writeNETCDFdata(n,ftn,CLM45_struc(n)%EF1_FDT,"double")
!      print*,'wrote EF1_FDT'
      call LDT_writeNETCDFdata(n,ftn,CLM45_struc(n)%EF1_SHR,"double")
!      print*,'wrote EF1_SHR'
      call LDT_writeNETCDFdata(n,ftn,CLM45_struc(n)%EF1_GRS,"double")
!      print*,'wrote EF1_GRS'
      call LDT_writeNETCDFdata(n,ftn,CLM45_struc(n)%EF1_CRP,"double")
!      print*,'wrote EF1_CRP'
      call LDT_writeNETCDFdata(n,ftn,CLM45_struc(n)%CANYON_HWR,"double")
!      print*,'wrote CANYON_HWR'
      call LDT_writeNETCDFdata(n,ftn,CLM45_struc(n)%EM_IMPROAD,"double")
!      print*,'wrote EM_IMPROAD'
      call LDT_writeNETCDFdata(n,ftn,CLM45_struc(n)%EM_PERROAD,"double")
!      print*,'wrote EM_PERROAD'
      call LDT_writeNETCDFdata(n,ftn,CLM45_struc(n)%EM_ROOF,"double")
!      print*,'wrote EM_ROOF'
      call LDT_writeNETCDFdata(n,ftn,CLM45_struc(n)%EM_WALL,"double")
!      print*,'wrote EM_WALL'
      call LDT_writeNETCDFdata(n,ftn,CLM45_struc(n)%HT_ROOF,"double")
!      print*,'wrote HT_ROOF'
      call LDT_writeNETCDFdata(n,ftn,CLM45_struc(n)%THICK_ROOF,"double")
!      print*,'wrote THICK_ROOF'
      call LDT_writeNETCDFdata(n,ftn,CLM45_struc(n)%THICK_WALL,"double")
!      print*,'wrote THICK_WALL'
      call LDT_writeNETCDFdata(n,ftn,CLM45_struc(n)%T_BUILDING_MAX,"double")
!      print*,'wrote T_BUILDING_MAX'
      call LDT_writeNETCDFdata(n,ftn,CLM45_struc(n)%T_BUILDING_MIN,"double")
!      print*,'wrote T_BUILDING_MIN'
      call LDT_writeNETCDFdata(n,ftn,CLM45_struc(n)%WIND_HGT_CANYON,"double")
!      print*,'wrote WIND_HGT_CANYON' 
      call LDT_writeNETCDFdata(n,ftn,CLM45_struc(n)%WTLUNIT_ROOF,"double")
!      print*,'wrote WTLUNIT_ROOF' 
      call LDT_writeNETCDFdata(n,ftn,CLM45_struc(n)%WTROAD_PERV,"double")
!      print*,'wrote WTROAD_PERV' 
      call LDT_writeNETCDFdata(n,ftn,CLM45_struc(n)%ALB_IMPROAD_DIR,"double")
!      print*,'wrote ALB_IMPROAD_DIR' 
      call LDT_writeNETCDFdata(n,ftn,CLM45_struc(n)%ALB_IMPROAD_DIF,"double")
!      print*,'wrote ALB_IMPROAD_DIF' 
      call LDT_writeNETCDFdata(n,ftn,CLM45_struc(n)%ALB_PERROAD_DIR,"double")
!      print*,'wrote ALB_PERROAD_DIR' 
      call LDT_writeNETCDFdata(n,ftn,CLM45_struc(n)%ALB_PERROAD_DIF,"double")
!      print*,'wrote ALB_PERROAD_DIF' 
      call LDT_writeNETCDFdata(n,ftn,CLM45_struc(n)%ALB_ROOF_DIR,"double")
!      print*,'wrote ALB_ROOF_DIR' 
      call LDT_writeNETCDFdata(n,ftn,CLM45_struc(n)%ALB_ROOF_DIF,"double")
!      print*,'wrote ALB_ROOF_DIF' 
      call LDT_writeNETCDFdata(n,ftn,CLM45_struc(n)%ALB_WALL_DIR,"double")
!      print*,'wrote ALB_WALL_DIR' 
      call LDT_writeNETCDFdata(n,ftn,CLM45_struc(n)%ALB_WALL_DIF,"double")
!      print*,'wrote ALB_WALL_DIF' 
      call LDT_writeNETCDFdata(n,ftn,CLM45_struc(n)%TK_ROOF,"double")
!      print*,'wrote TK_ROOF' 
      call LDT_writeNETCDFdata(n,ftn,CLM45_struc(n)%TK_WALL,"double")
!      print*,'wrote TK_WALL' 
      call LDT_writeNETCDFdata(n,ftn,CLM45_struc(n)%TK_IMPROAD,"double")
!      print*,'wrote TK_IMPROAD' 
      call LDT_writeNETCDFdata(n,ftn,CLM45_struc(n)%CV_ROOF,"double")
!      print*,'wrote CV_ROOF' 
      call LDT_writeNETCDFdata(n,ftn,CLM45_struc(n)%CV_WALL,"double")
!      print*,'wrote CV_WALL' 
      call LDT_writeNETCDFdata(n,ftn,CLM45_struc(n)%CV_IMPROAD,"double")
!      print*,'wrote CV_IMPROAD' 
      call LDT_writeNETCDFdata(n,ftn,CLM45_struc(n)%NLEV_IMPROAD,"int")
!      print*,'wrote NLEV_IMPROAD' 
      call LDT_writeNETCDFdata(n,ftn,CLM45_struc(n)%peatf,"double")
!      print*,'wrote peatf' 
      call LDT_writeNETCDFdata(n,ftn,CLM45_struc(n)%abm,"int")
!      print*,'wrote abm' 
      call LDT_writeNETCDFdata(n,ftn,CLM45_struc(n)%gdp,"double")
!      print*,'wrote gdp' 
      call LDT_writeNETCDFdata(n,ftn,CLM45_struc(n)%SLOPE,"double")
!      print*,'wrote SLOPE' 
      call LDT_writeNETCDFdata(n,ftn,CLM45_struc(n)%STD_ELEV,"double")
!      print*,'wrote STD_ELEV' 
      call LDT_writeNETCDFdata(n,ftn,CLM45_struc(n)%binfl,"double")
!      print*,'wrote binfl' 
      call LDT_writeNETCDFdata(n,ftn,CLM45_struc(n)%Ws,"double")
!      print*,'wrote Ws' 
      call LDT_writeNETCDFdata(n,ftn,CLM45_struc(n)%Dsmax,"double")
!      print*,'wrote Dsmax' 
      call LDT_writeNETCDFdata(n,ftn,CLM45_struc(n)%Ds,"double")
!      print*,'wrote Ds' 
      call LDT_writeNETCDFdata(n,ftn,CLM45_struc(n)%LAKEDEPTH,"double")
!      print*,'wrote LAKEDEPTH' 
      call LDT_writeNETCDFdata(n,ftn,CLM45_struc(n)%F0,"double")
!      print*,'wrote F0' 
      call LDT_writeNETCDFdata(n,ftn,CLM45_struc(n)%P3,"double")
!      print*,'wrote P3' 
      call LDT_writeNETCDFdata(n,ftn,CLM45_struc(n)%ZWT0,"double")
!      print*,'wrote ZWT0' 
      call LDT_writeNETCDFdata(n,ftn,CLM45_struc(n)%PCT_WETLAND,"double")
!      print*,'wrote PCT_WETLAND' 
      call LDT_writeNETCDFdata(n,ftn,CLM45_struc(n)%PCT_LAKE,"double")
!      print*,'wrote PCT_LAKE' 
      call LDT_writeNETCDFdata(n,ftn,CLM45_struc(n)%PCT_GLACIER,"double")
!      print*,'wrote PCT_GLACIER' 
      ! special 1D variable: GLC_MEC
      call LDT_verify(nf90_put_var(ftn, CLM45_struc(n)%GLC_MEC%vid,&
           CLM45_struc(n)%GLC_MEC%dvalue(:,1,1,1), &
           (/1,1,1,1/),(/CLM45_struc(n)%nglcecp1,1,1,1/)), &
           'nf90_put_var failed in CLM45parms_writeData')
         call LDT_verify(nf90_put_att(ftn,CLM45_struc(n)%GLC_MEC%vid,&
              "standard_name",trim(CLM45_struc(n)%GLC_MEC%standard_name)))
         call LDT_verify(nf90_put_att(ftn,CLM45_struc(n)%GLC_MEC%vid,&
              "units",trim(CLM45_struc(n)%GLC_MEC%units)))
      call LDT_writeNETCDFdata(n,ftn,CLM45_struc(n)%PCT_GLC_MEC,"double")
!      print*,'wrote PCT_GLC_MEC' 
      call LDT_writeNETCDFdata(n,ftn,CLM45_struc(n)%PCT_GLC_MEC_GIC,"double")
!      print*,'wrote PCT_GLC_MEC_GIC' 
      call LDT_writeNETCDFdata(n,ftn,CLM45_struc(n)%PCT_GLC_MEC_ICESHEET,"double")
!      print*,'wrote PCT_GLC_MEC_ICESHEET' 
      call LDT_writeNETCDFdata(n,ftn,CLM45_struc(n)%PCT_GLC_GIC,"double")
!      print*,'wrote PCT_GLC_GIC' 
      call LDT_writeNETCDFdata(n,ftn,CLM45_struc(n)%PCT_GLC_ICESHEET,"double")
!      print*,'wrote PCT_GLC_ICESHEET' 
      call LDT_writeNETCDFdata(n,ftn,CLM45_struc(n)%TOPO_GLC_MEC,"double")
!      print*,'wrote TOPO_GLC_MEC' 
      call LDT_writeNETCDFdata(n,ftn,CLM45_struc(n)%TOPO,"double")
!      print*,'wrote TOPO' 
      call LDT_writeNETCDFdata(n,ftn,CLM45_struc(n)%PCT_URBAN,"double")
!      print*,'wrote PCT_URBAN' 
      call LDT_writeNETCDFdata(n,ftn,CLM45_struc(n)%URBAN_REGION_ID,"int")
!      print*,'wrote URBAN_REGION_ID' 

    endif

  end subroutine CLM45parms_writeData

!BOP
! !ROUTINE:  set_param_attribs
! \label{set_param_attribs}
!
! !INTERFACE:
  subroutine set_param_attribs(paramEntry, short_name, long_name, units, &
                               vlevels, zlevels)

! !DESCRIPTION:
!   This routine reads over the parameter attribute entries
!   in the param_attribs.txt file.
!
! !USES:
   type(LDT_paramEntry),intent(inout) :: paramEntry
   character(len=*),    intent(in)    :: short_name
   character(len=*),    intent(in)    :: long_name
   character(len=*),    intent(in)    :: units
   integer, optional                  :: vlevels
   integer, optional                  :: zlevels

   integer   :: v_temp, z_temp
! ____________________________________________________
    
   if(present(vlevels)) then
      v_temp = vlevels
   else
      v_temp = 1
   endif
   if(present(zlevels)) then
      z_temp = zlevels
   else
      z_temp = 1
   endif

! LIS starndard addition
!   paramEntry%long_name =trim(long_name)
!   paramEntry%scale_factor = 1
!   paramEntry%add_offset = 0.
   paramEntry%short_name = trim(short_name)
   paramEntry%vlevels = v_temp
   paramEntry%zlevels = z_temp
   paramEntry%selectOpt = 1
   paramEntry%source = "CLM.4.5"
   paramEntry%num_times = 1
   paramEntry%num_bins = 1
   paramEntry%units =trim(units)
   paramEntry%standard_name = trim(long_name)

  end subroutine set_param_attribs

end module CLM45_parmsMod

