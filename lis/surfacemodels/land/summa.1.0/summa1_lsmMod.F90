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
module summa1_lsmMod
!BOP
!
! !MODULE: summa1_lsmMod
!
! !DESCRIPTION:
!  Module for 1-D land model driver variable initialization
!  
! \begin{description}
!  \item[count]
!    variable to keep track of the number of timesteps before an output
!  \item[numout]
!    number of output times 
!  \item[outInterval]
!    output writing interval
!  \item[summa1open]
!    variable to keep track of opened files
!  \item[summa1]
!   Summa1 LSM specific variables
! \end{description} 
!
! !REVISION HISTORY:
!   1 Jun 2016; Sujay Kumar, Initial Code
!
! !USES:        
  use nrtype
  use summa1_module
  use globalData

  implicit none
  
  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: summa1_lsm_ini
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  public :: summa1_struc
!EOP
  type, public :: summa1_type_dec
! *****************************************************************************
! (0) variable definitions
! *****************************************************************************
     type(gru_hru_doubleVec)          :: forcStat                   ! x%gru(:)%hru(:)%var(:)%dat -- model forcing data
     type(gru_hru_doubleVec)          :: progStat                   ! x%gru(:)%hru(:)%var(:)%dat -- model prognostic (state) variables
     type(gru_hru_doubleVec)          :: diagStat                   ! x%gru(:)%hru(:)%var(:)%dat -- model diagnostic variables
     type(gru_hru_doubleVec)          :: fluxStat                   ! x%gru(:)%hru(:)%var(:)%dat -- model fluxes
     type(gru_hru_doubleVec)          :: indxStat                   ! x%gru(:)%hru(:)%var(:)%dat -- model indices
     type(gru_doubleVec)              :: bvarStat                   ! x%gru(:)%var(:)%dat        -- basin-average variabl
     ! define the primary data structures (scalars)
     type(var_i)                      :: timeStruct                 ! x%var(:)                   -- model time data
     type(gru_hru_double)             :: forcStruct                 ! x%gru(:)%hru(:)%var(:)     -- model forcing data
     type(gru_hru_double)             :: attrStruct                 ! x%gru(:)%hru(:)%var(:)     -- local attributes for each HRU
     type(gru_hru_int)                :: typeStruct                 ! x%gru(:)%hru(:)%var(:)     -- local classification of soil veg etc. for each HRU
     type(gru_hru_doubleVec)          :: mparStruct                 ! x%gru(:)%hru(:)%var(:)     -- model parameters
     ! define the primary data structures (variable length vectors)
     type(gru_hru_intVec)             :: indxStruct                 ! x%gru(:)%hru(:)%var(:)%dat -- model indices
     type(gru_hru_doubleVec)          :: progStruct                 ! x%gru(:)%hru(:)%var(:)%dat -- model prognostic (state) variables
     type(gru_hru_doubleVec)          :: diagStruct                 ! x%gru(:)%hru(:)%var(:)%dat -- model diagnostic variables
     type(gru_hru_doubleVec)          :: fluxStruct                 ! x%gru(:)%hru(:)%var(:)%dat -- model fluxes
     ! define the basin-average structures
     type(gru_double)                 :: bparStruct                 ! x%gru(:)%var(:)            -- basin-average parameters
     type(gru_doubleVec)              :: bvarStruct                 ! x%gru(:)%var(:)%dat        -- basin-average variables
     ! define the ancillary data structures
     type(gru_hru_double)             :: dparStruct                 ! x%gru(:)%hru(:)%var(:)     -- default model parameters
     
     type(hru_i),allocatable          :: computeVegFlux(:)          ! flag to indicate if we are computing fluxes over vegetation (.false. means veg is buried with snow)      
     type(hru_d),allocatable          :: dt_init(:)                 ! used to initialize the length of the sub-step for each HRU
     type(hru_d),allocatable          :: upArea(:)                  ! area upslope of each HRU 

     integer                    :: summa1open
     integer                    :: numout
     real                       :: ts
     integer(i4b)                     :: nGRU                       ! number of grouped response units
     integer(i4b)                     :: nHRU                       ! number of global hydrologic response units
     integer(i4b)                     :: hruCount                   ! number of local hydrologic response units
     real(dp),dimension(12)           :: greenVegFrac_monthly       ! fraction of green vegetation in each month (0-1)
     character(len=256)               :: summaFileManagerFile   ! path/name of file defining directories and files

     type(summa1dec), allocatable :: summa1(:)
  end type summa1_type_dec

  type(summa1_type_dec), allocatable :: summa1_struc(:)

  SAVE
contains
!BOP
! 
! !ROUTINE: summa1_lsm_ini
! \label{summa1_lsm_ini}
! 
! !INTERFACE:
  subroutine summa1_lsm_ini()
! !USES:
   use ESMF
   use LIS_surfaceModelDataMod, only : LIS_sfmodel_struc
   use LIS_coreMod, only : LIS_rc
   use LIS_timeMgrMod,   only : LIS_clock, LIS_calendar, &
        LIS_update_timestep, LIS_registerAlarm
   use LIS_logMod,       only : LIS_verify
   
   USE allocspace_module
   USE summaFileManager
   USE popMetadat_module
   USE checkStruc_module
   USE childStruc_module
   USE var_lookup
   USE ascii_util_module
   USE read_attrb_module
   USE ffile_info_module
   USE mDecisions_module
   USE read_pinit_module
   USE module_sf_noahmplsm
   USE pOverwrite_module
   USE read_param_module
   USE var_derive_module
   USE paramCheck_module
   USE convE2Temp_module
   USE read_icond_module
   USE flxMapping_module
   USE output_stats
   USE modelwrite_module
   USE def_output_module
   USE check_icond_module
! !DESCRIPTION:        
!
!EOP
   implicit none
   integer :: n
   integer                 :: yr, mo, da, hr, mn, ss
   integer                 :: status

   integer(i4b)                     :: err=0                      ! error code
   character(len=1024)              :: message=''                 ! error message
!   character(len=256)               :: summaFileManagerFile=''    ! path/name of file defining directories and files
   logical(lgt)                     :: flux_mask(size(flux_meta)) ! mask defining desired flux variables
!   integer(i4b)                     :: nGRU                       ! number of grouped response units
!   integer(i4b)                     :: nHRU                       ! number of global hydrologic response units
   integer(i4b)                     :: hruCount                   ! number of local hydrologic response units
   integer(i4b)                     :: iVar                       ! index of a model variable 
   integer(i4b)                     :: iGRU
   integer(i4b)                     :: iHRU,jHRU,kHRU             ! index of the hydrologic response unit
   integer(i4b)                     :: iWord                      ! loop through words in a string
   integer(i4b)                     :: iStruct                    ! loop through data structures
   integer(i4b)                     :: fileUnit                   ! file unit (output from file_open; a unit not currently used)
   character(LEN=256),allocatable   :: dataLines(:)               ! vector of character strings from non-comment lines
   character(LEN=256),allocatable   :: chardata(:)                ! vector of character data
   character(len=256)               :: fileout=''                 ! output filename
   character(len=64)                :: output_fileSuffix=''       ! suffix for the output file
   character(len=256)               :: attrFile                   ! attributes file name
   character(len=256)               :: restartFile                ! restart file name
   integer(i4b)                     :: iRunMode                   ! define the current running mode
   integer(i4b)                     :: fileGRU                    ! number of GRUs in the input file
   integer(i4b)                     :: fileHRU                    ! number of HRUs in the input file
   integer(i4b)                     :: startGRU                   ! index of the starting GRU for parallelization run
   integer(i4b)                     :: checkHRU                   ! index of the HRU for a single HRU run

   allocate(summa1_struc(LIS_rc%nnest))

   call summa1_readcrd()
   do n = 1, LIS_rc%nnest
      allocate(summa1_struc(n)%summa1(LIS_rc%npatch(n,LIS_rc%lsm_index)))
      summa1_struc(n)%numout = 0

      call LIS_update_timestep(LIS_rc, n, summa1_struc(n)%ts)

      LIS_sfmodel_struc(n)%nsm_layers = 1
      LIS_sfmodel_struc(n)%nst_layers = 1
      allocate(LIS_sfmodel_struc(n)%lyrthk(1))
      LIS_sfmodel_struc(n)%lyrthk(1) = 1
      LIS_sfmodel_struc(n)%ts = summa1_struc(n)%ts
   enddo

   do n = 1, LIS_rc%nnest

!initializing variables, from getCommandArguments
      startGRU = integerMissing; 
      checkHRU = integerMissing; 
      summa1_struc(n)%nHRU = integerMissing
      summa1_struc(n)%nGRU = integerMissing; 
      iRunMode = iRunModeFull
      
!TBD - hardcoded 
!   summaFileManagerFile = './summaTestCases/settings_mergeAlloc/syntheticTestCases/celia1990/summa_fileManager_celia1990.txt'
      call summa_SetDirsUndPhiles(summa1_struc(n)%summaFileManagerFile,err,message)
      
      call allocLocal(time_meta, refTime,   err=err, message=message)
      call allocLocal(time_meta, startTime, err=err, message=message)
      call allocLocal(time_meta, finshTime, err=err, message=message)

! *****************************************************************************
! (2) populate/check metadata structures
! *****************************************************************************  
! populate metadata for all model variables
      call popMetadat(err,message);    call LIS_verify(err,message)
! define mapping between fluxes and states
      call flxMapping(err,message);   call LIS_verify(err,message)
! check data structures
      call checkStruc(err,message);   call LIS_verify(err,message)

! define the mask to identify the subset of variables in the "child" data structure (just scalar variables)
      flux_mask = (flux_meta(:)%vartype==iLookVarType%scalarv)

   ! create the averageFlux metadata structure
      call childStruc(flux_meta, flux_mask, averageFlux_meta, childFLUX_MEAN, err, message);   call LIS_verify(err,message)

! *****************************************************************************
! (3a) read the number of GRUs and HRUs, and allocate the gru-hru mapping structures
! *****************************************************************************

      attrFile = trim(SETNGS_PATH)//trim(LOCAL_ATTRIBUTES)
      
      print*, 'here',attrFile, iRunMode, iRunMode, iRunModeGRU, iRunModeHRU
      fileGRU = LIS_rc%ngrid(n)
      fileHRU = LIS_rc%ngrid(n)
      summa1_struc(n)%nGRU = LIS_rc%ngrid(n)
      summa1_struc(n)%nHRU = LIS_rc%ngrid(n)
      call read_dimension_LIS(trim(LIS_rc%paramfile(n)),fileGRU,fileHRU,&
           summa1_struc(n)%nGRU,summa1_struc(n)%nHRU,err,message)
! Commented out by SVK     
!      select case (iRunMode)
!      case(iRunModeFull); call read_dimension(trim(LIS_rc%paramfile),fileGRU,fileHRU,&
!           summa1_struc(n)%nGRU,summa1_struc(n)%nHRU,err,message)
!      case(iRunModeGRU ); call read_dimension(trim(attrFile),fileGRU,fileHRU,&
!           summa1_struc(n)%nGRU,summa1_struc(n)%nHRU,err,message,startGRU=startGRU)
!      case(iRunModeHRU ); call read_dimension(trim(attrFile),fileGRU,fileHRU,&
!           summa1_struc(n)%nGRU,summa1_struc(n)%nHRU,err,message,checkHRU=checkHRU)
!      end select
      call LIS_verify(err,message)
      print*, fileGRU, fileHRU, summa1_struc(n)%nGRU, summa1_struc(n)%nHRU

! *****************************************************************************
! (3b) read model attributes 
! *****************************************************************************
! SVK: This file reading needs to be modified for restarts - Fix it 
! when proper restart mode is implemented. Currently this will not 
! work for parallel runs. 

! read number of snow and soil layers
      restartFile = trim(SETNGS_PATH)//trim(MODEL_INITCOND)
      print*, trim(restartFile)
      call read_icond_nlayers(trim(restartFile),summa1_struc(n)%nGRU,indx_meta,err,message)
      call LIS_verify(err,message)

! *****************************************************************************
! (3c) allocate space for other data structures
! *****************************************************************************
! loop through data structures
   do iStruct=1,size(structInfo)
      ! allocate space
      select case(trim(structInfo(iStruct)%structName))
      case('time'); call allocGlobal(time_meta,  summa1_struc(n)%timeStruct,  err, message)   ! model forcing data
      case('forc'); call allocGlobal(forc_meta,  summa1_struc(n)%forcStruct,  err, message)   ! model forcing data
      case('attr'); call allocGlobal(attr_meta,  summa1_struc(n)%attrStruct,  err, message)   ! local attributes for each HRU
      case('type'); call allocGlobal(type_meta,  summa1_struc(n)%typeStruct,  err, message)   ! local classification of soil veg etc. for each HRU
      case('mpar'); call allocGlobal(mpar_meta,  summa1_struc(n)%mparStruct,  err, message)   ! model parameters
      case('indx'); call allocGlobal(indx_meta,  summa1_struc(n)%indxStruct,  err, message)   ! model variables
      case('prog'); call allocGlobal(prog_meta,  summa1_struc(n)%progStruct,  err, message)   ! model prognostic (state) variables
      case('diag'); call allocGlobal(diag_meta,  summa1_struc(n)%diagStruct,  err, message)   ! model diagnostic variables
      case('flux'); call allocGlobal(flux_meta,  summa1_struc(n)%fluxStruct,  err, message)   ! model fluxes
      case('bpar'); call allocGlobal(bpar_meta,  summa1_struc(n)%bparStruct,  err, message)   ! basin-average parameters
      case('bvar'); call allocGlobal(bvar_meta,  summa1_struc(n)%bvarStruct,  err, message)   ! basin-average variables
      case('deriv'); cycle
      case default; err=20; message='unable to find structure name: '//trim(structInfo(iStruct)%structName)
      end select
      ! check errors
      call LIS_verify(err,trim(message)//'[structure =  '//trim(structInfo(iStruct)%structName)//']')
   end do  ! looping through data structures
   
! *****************************************************************************
! (3c) allocate space for other data structures
! allocate space for default model parameters
! NOTE: This is done here, rather than in the loop above, because dpar is not one of the "standard" data structures
! *****************************************************************************
   call allocGlobal(mpar_meta,summa1_struc(n)%dparStruct,err,message)   ! default model parameters
   call LIS_verify(err,trim(message)//' [problem allocating dparStruct]')
   
   allocate(summa1_struc(n)%dt_init(summa1_struc(n)%nGRU),&
        summa1_struc(n)%upArea(summa1_struc(n)%nGRU),&
        summa1_struc(n)%computeVegFlux(summa1_struc(n)%nGRU),stat=err)
   call LIS_verify(err,'problem allocating space for dt_init, upArea, or computeVegFlux [GRU]')

! allocate space for the HRUs
   do iGRU=1,summa1_struc(n)%nGRU
      hruCount = gru_struc(iGRU)%hruCount
      allocate(summa1_struc(n)%dt_init(iGRU)%hru(hruCount),&
           summa1_struc(n)%upArea(iGRU)%hru(hruCount),&
           summa1_struc(n)%computeVegFlux(iGRU)%hru(hruCount),stat=err)
      call LIS_verify(err,'problem allocating space for dt_init, upArea, or computeVegFlux [HRU]')
   end do

! *****************************************************************************
! (4a) read local attributes for each HRU
! *****************************************************************************
   call read_attrb_LIS(n,trim(LIS_rc%paramfile(n)),&
        summa1_struc(n)%nGRU,&
        summa1_struc(n)%attrStruct,&
        summa1_struc(n)%typeStruct,err,message)
   call LIS_verify(err,message)

! *****************************************************************************
! (4b) read description of model forcing datafile used in each HRU
! *****************************************************************************
!call ffile_info(nGRU,err,message); call LIS_verify(err,message)
   data_step = summa1_struc(n)%ts
! *****************************************************************************
! (4c) read model decisions
! *****************************************************************************
   call mDecisions(err,message); call LIS_verify(err,message)


! *****************************************************************************
! (5a) read default model parameters
! *****************************************************************************
! read default values and constraints for model parameters (local column, and basin-average)
   call read_pinit(LOCALPARAM_INFO,.TRUE., mpar_meta,localParFallback,err,message)
   call LIS_verify(err,message)
   call read_pinit(BASINPARAM_INFO,.FALSE.,bpar_meta,basinParFallback,err,message); 
   call LIS_verify(err,message)


! *****************************************************************************
! (5b) read Noah vegetation and soil tables
! *****************************************************************************
! define monthly fraction of green vegetation
   summa1_struc(n)%greenVegFrac_monthly = (/0.01_dp, 0.02_dp, 0.03_dp, 0.07_dp, 0.50_dp, 0.90_dp, 0.95_dp, 0.96_dp, 0.65_dp, 0.24_dp, 0.11_dp, 0.02_dp/)
   
   ! read Noah soil and vegetation tables
   call soil_veg_gen_parm(trim(SETNGS_PATH)//'VEGPARM.TBL',                              & ! filename for vegetation table
        trim(SETNGS_PATH)//'SOILPARM.TBL',                             & ! filename for soils table
        trim(SETNGS_PATH)//'GENPARM.TBL',                              & ! filename for general table
        trim(model_decisions(iLookDECISIONS%vegeParTbl)%cDecision),    & ! classification system used for vegetation
        trim(model_decisions(iLookDECISIONS%soilCatTbl)%cDecision))      ! classification system used for soils
   
! read Noah-MP vegetation tables
   call read_mp_veg_parameters(trim(SETNGS_PATH)//'MPTABLE.TBL',                         & ! filename for Noah-MP table
        trim(model_decisions(iLookDECISIONS%vegeParTbl)%cDecision)) ! classification system used for vegetation
   
! define urban vegetation category
   select case(trim(model_decisions(iLookDECISIONS%vegeParTbl)%cDecision))
   case('USGS');                     urbanVegCategory =    1
   case('MODIFIED_IGBP_MODIS_NOAH'); urbanVegCategory =   13
   case('plumberCABLE');             urbanVegCategory = -999
   case('plumberCHTESSEL');          urbanVegCategory = -999
   case('plumberSUMMA');             urbanVegCategory = -999
   case default; call LIS_verify(30,'unable to identify vegetation category')
   end select

! set default model parameters
   do iGRU=1, summa1_struc(n)%nGRU
      do iHRU=1,gru_struc(iGRU)%hruCount
         ! set parmameters to their default value
          summa1_struc(n)%dparStruct%gru(iGRU)%hru(iHRU)%var(:) = localParFallback(:)%default_val         ! x%hru(:)%var(:)
         ! overwrite default model parameters with information from the Noah-MP tables
         call pOverwrite( summa1_struc(n)%typeStruct%gru(iGRU)%hru(iHRU)%var(iLookTYPE%vegTypeIndex),  &  ! vegetation category
               summa1_struc(n)%typeStruct%gru(iGRU)%hru(iHRU)%var(iLookTYPE%soilTypeIndex), &  ! soil category
               summa1_struc(n)%dparStruct%gru(iGRU)%hru(iHRU)%var,                          &  ! default model parameters
              err,message); call LIS_verify(err,message)            ! error control
         ! copy over to the parameter structure
         ! NOTE: constant for the dat(:) dimension (normally depth)
         do ivar=1,size(localParFallback)
             summa1_struc(n)%mparStruct%gru(iGRU)%hru(iHRU)%var(ivar)%dat(:) = &
                  summa1_struc(n)%dparStruct%gru(iGRU)%hru(iHRU)%var(ivar)
         end do  ! looping through variables
      end do  ! looping through HRUs
      ! set default for basin-average parameters
      summa1_struc(n)%bparStruct%gru(iGRU)%var(:) = basinParFallback(:)%default_val
   end do  ! looping through GRUs

! *****************************************************************************
! (5c) read trial model parameter values for each HRU, and populate initial data structures
! *****************************************************************************
   call read_param(iRunMode,checkHRU,startGRU, summa1_struc(n)%nHRU, summa1_struc(n)%nGRU,&
        summa1_struc(n)%typeStruct,&
        summa1_struc(n)%mparStruct, &
        summa1_struc(n)%bparStruct,&
        err,message); 
   call LIS_verify(err,message)
   

! *****************************************************************************
! (5d) compute derived model variables that are pretty much constant for the basin as a whole
! *****************************************************************************
! loop through GRUs
   do iGRU=1,summa1_struc(n)%nGRU

 ! calculate the fraction of runoff in future time steps
      call fracFuture(summa1_struc(n)%bparStruct%gru(iGRU)%var,    &  ! vector of basin-average model parameters
           summa1_struc(n)%bvarStruct%gru(iGRU),        &  ! data structure of basin-average variables
           err,message)                    ! error control
      call LIS_verify(err,message)
      
 ! loop through local HRUs
      do iHRU=1,gru_struc(iGRU)%hruCount
         
         kHRU=0
         ! check the network topology (only expect there to be one downslope HRU)
         do jHRU=1,gru_struc(iGRU)%hruCount
            if(summa1_struc(n)%typeStruct%gru(iGRU)%hru(iHRU)%var(iLookTYPE%downHRUindex) == &
                 summa1_struc(n)%typeStruct%gru(iGRU)%hru(jHRU)%var(iLookTYPE%hruIndex))then
               if(kHRU==0)then  ! check there is a unique match
                  kHRU=jHRU
               else
                  call LIS_verify(20,'multi_driver: only expect there to be one downslope HRU')
               end if  ! (check there is a unique match)
            end if  ! (if identified a downslope HRU)
         end do
  
  ! check that the parameters are consistent
         call paramCheck(summa1_struc(n)%mparStruct%gru(iGRU)%hru(iHRU),err,message); call LIS_verify(err,message)
  
  ! calculate a look-up table for the temperature-enthalpy conversion 
         call E2T_lookup(summa1_struc(n)%mparStruct%gru(iGRU)%hru(iHRU),err,message); call LIS_verify(err,message)
         
      end do ! HRU
   end do ! GRU


! read description of model initial conditions -- also initializes model structure components
! NOTE: at this stage the same initial conditions are used for all HRUs -- need to modify
   call read_icond(restartFile,                   & ! name of initial conditions file
        summa1_struc(n)%nGRU,                          & ! number of response units
        prog_meta,                     & ! metadata
        summa1_struc(n)%progStruct,                    & ! model prognostic (state) variables
        summa1_struc(n)%indxStruct,                    & ! layer indexes
        err,message)                     ! error control
   call LIS_verify(err,message)
   
! check initial conditions
   call check_icond(summa1_struc(n)%nGRU,                          & ! number of response units
        summa1_struc(n)%progStruct,                    & ! model prognostic (state) variables
        summa1_struc(n)%mparStruct,                    & ! model parameters
        summa1_struc(n)%indxStruct,                    & ! layer indexes
        err,message)                     ! error control
   call LIS_verify(err,message)

! loop through GRUs
   do iGRU=1, summa1_struc(n)%nGRU
 ! loop through local HRUs
      do iHRU=1,gru_struc(iGRU)%hruCount

  ! re-calculate height of each layer
         call calcHeight(&
              ! input/output: data structures
              summa1_struc(n)%indxStruct%gru(iGRU)%hru(iHRU),   & ! intent(in): layer type
              summa1_struc(n)%progStruct%gru(iGRU)%hru(iHRU),   & ! intent(inout): model prognostic (state) variables for a local HRU
              ! output: error control
              err,message); call LIS_verify(err,message)
         
         ! calculate vertical distribution of root density
         call rootDensty( summa1_struc(n)%mparStruct%gru(iGRU)%hru(iHRU),    & ! vector of model parameters
               summa1_struc(n)%indxStruct%gru(iGRU)%hru(iHRU),    & ! data structure of model indices
               summa1_struc(n)%progStruct%gru(iGRU)%hru(iHRU),    & ! data structure of model prognostic (state) variables
               summa1_struc(n)%diagStruct%gru(iGRU)%hru(iHRU),    & ! data structure of model diagnostic variables
              err,message)                         ! error control
         call LIS_verify(err,message) 
         
         ! calculate saturated hydraulic conductivity in each soil layer
         call satHydCond( summa1_struc(n)%mparStruct%gru(iGRU)%hru(iHRU),    & ! vector of model parameters
               summa1_struc(n)%indxStruct%gru(iGRU)%hru(iHRU),    & ! data structure of model indices
               summa1_struc(n)%progStruct%gru(iGRU)%hru(iHRU),    & ! data structure of model prognostic (state) variables
               summa1_struc(n)%fluxStruct%gru(iGRU)%hru(iHRU),    & ! data structure of model fluxes 
              err,message)                         ! error control
         call LIS_verify(err,message)
         
         ! calculate "short-cut" variables such as volumetric heat capacity
         call v_shortcut( summa1_struc(n)%mparStruct%gru(iGRU)%hru(iHRU),    & ! vector of model parameters
               summa1_struc(n)%diagStruct%gru(iGRU)%hru(iHRU),    & ! data structure of model diagnostic variables
              err,message)                         ! error control
         call LIS_verify(err,message)
         
         ! overwrite the vegetation height
         HVT(summa1_struc(n)%typeStruct%gru(iGRU)%hru(iHRU)%var(iLookTYPE%vegTypeIndex)) = &
              summa1_struc(n)%mparStruct%gru(iGRU)%hru(iHRU)%var(iLookPARAM%heightCanopyTop)%dat(1)
         HVB(summa1_struc(n)%typeStruct%gru(iGRU)%hru(iHRU)%var(iLookTYPE%vegTypeIndex)) = &
              summa1_struc(n)%mparStruct%gru(iGRU)%hru(iHRU)%var(iLookPARAM%heightCanopyBottom)%dat(1)
         
         ! overwrite the tables for LAI and SAI
         if(model_decisions(iLookDECISIONS%LAI_method)%iDecision == specified)then
            SAIM( summa1_struc(n)%typeStruct%gru(iGRU)%hru(iHRU)%var(iLookTYPE%vegTypeIndex),:) =&
                 summa1_struc(n)%mparStruct%gru(iGRU)%hru(iHRU)%var(iLookPARAM%winterSAI)%dat(1)
            LAIM( summa1_struc(n)%typeStruct%gru(iGRU)%hru(iHRU)%var(iLookTYPE%vegTypeIndex),:) = &
                 summa1_struc(n)%mparStruct%gru(iGRU)%hru(iHRU)%var(iLookPARAM%summerLAI)%dat(1)*&
                 summa1_struc(n)%greenVegFrac_monthly
         endif
         
         ! initialize canopy drip
  ! NOTE: canopy drip from the previous time step is used to compute throughfall for the current time step
          summa1_struc(n)%fluxStruct%gru(iGRU)%hru(iHRU)%var(iLookFLUX%scalarCanopyLiqDrainage)%dat(1) = 0._dp  ! not used  
      end do  ! (looping through HRUs)
      
      ! compute total area of the upstream HRUS that flow into each HRU
      do iHRU=1,gru_struc(iGRU)%hruCount
         summa1_struc(n)%upArea(iGRU)%hru(iHRU) = 0._dp
         do jHRU=1,gru_struc(iGRU)%hruCount
            ! check if jHRU flows into iHRU; assume no exchange between GRUs
            if( summa1_struc(n)%typeStruct%gru(iGRU)%hru(jHRU)%var(iLookTYPE%downHRUindex)==&
                  summa1_struc(n)%typeStruct%gru(iGRU)%hru(iHRU)%var(iLookTYPE%hruIndex))then
               summa1_struc(n)%upArea(iGRU)%hru(iHRU) = summa1_struc(n)%upArea(iGRU)%hru(iHRU) +  &
                    summa1_struc(n)%attrStruct%gru(iGRU)%hru(jHRU)%var(iLookATTR%HRUarea)
            endif   ! (if jHRU is an upstream HRU)
         end do  ! jHRU
      end do  ! iHRU
      
      ! identify the total basin area for a GRU (m2)
      associate(totalArea => summa1_struc(n)%bvarStruct%gru(iGRU)%var(iLookBVAR%basin__totalArea)%dat(1) )
      totalArea = 0._dp
      do iHRU=1,gru_struc(iGRU)%hruCount
         totalArea = totalArea + summa1_struc(n)%attrStruct%gru(iGRU)%hru(iHRU)%var(iLookATTR%HRUarea)
      end do
      end associate

 ! initialize aquifer storage
 ! NOTE: this is ugly: need to add capabilities to initialize basin-wide state variables
 ! There are two options for groundwater:
 !  (1) where groundwater is included in the local column (i.e., the HRUs); and
 !  (2) where groundwater is included for the single basin (i.e., the GRUS, where multiple HRUS drain into a GRU).
 ! For water balance calculations it is important to ensure that the local aquifer storage is zero if groundwater is treated as a basin-average state variable (singleBasin);
 !  and ensure that basin-average aquifer storage is zero when groundwater is included in the local columns (localColumn).
      select case(model_decisions(iLookDECISIONS%spatial_gw)%iDecision)
  ! the basin-average aquifer storage is not used if the groundwater is included in the local column
      case(localColumn)
         summa1_struc(n)%bvarStruct%gru(iGRU)%var(iLookBVAR%basin__AquiferStorage)%dat(1) = 0._dp ! set to zero to be clear that there is no basin-average aquifer storage in this configuration 
  ! NOTE: the local column aquifer storage is not used if the groundwater is basin-average
  ! (i.e., where multiple HRUs drain to a basin-average aquifer)
      case(singleBasin)
         summa1_struc(n)%bvarStruct%gru(iGRU)%var(iLookBVAR%basin__AquiferStorage)%dat(1) = 1._dp
         do iHRU=1,gru_struc(iGRU)%hruCount
            summa1_struc(n)%progStruct%gru(iGRU)%hru(iHRU)%var(iLookPROG%scalarAquiferStorage)%dat(1) = 0._dp  ! set to zero to be clear that there is no local aquifer storage in this configuration
         end do
      case default; call LIS_verify(20,'unable to identify decision for regional representation of groundwater')
      end select

 ! initialize time step length for each HRU
      do iHRU=1,gru_struc(iGRU)%hruCount
         summa1_struc(n)%dt_init(iGRU)%hru(iHRU) =&
              summa1_struc(n)%progStruct%gru(iGRU)%hru(iHRU)%var(iLookPROG%dt_init)%dat(1) ! seconds
      end do

   end do  ! (looping through GRUs)

#if 0    
! *** TEMPORARY CODE ***
! code will be replaced once merge with the NetCDF branch

! NOTE: currently the same initial conditions for all HRUs and GRUs; will change when shift to NetCDF

! loop through GRUs
   do iGRU=1,summa1_struc(n)%nGRU

 ! check the GRU-HRU mapping structure is allocated 
      if(.not.allocated(gru_struc(iGRU)%hruInfo)) print*, 'gru_struc(iGRU)%hruInfo is not allocated'

 ! loop through HRUs
      do iHRU=1,gru_struc(iGRU)%hruCount  ! loop through HRUs within a given GRU

  ! get a vector of non-commented lines
         call file_open(trim(SETNGS_PATH)//trim(MODEL_INITCOND),fileUnit,err,message)
         call get_vlines(fileUnit,dataLines,err,message)
         close(fileUnit)
  
  ! get the number of snow and soil layers for each HRU
         gru_struc(iGRU)%hruInfo(iHRU)%nSnow = 0 ! initialize the number of snow layers
         gru_struc(iGRU)%hruInfo(iHRU)%nSoil = 0 ! initialize the number of soil layers
         do iVar=1,size(dataLines)
   ! split the line into an array of words
            call split_line(dataLines(iVar),chardata,err,message)
   ! check if the line contains initial conditions data (contains the word "snow" or "soil")
            do iword=1,size(chardata)
               if(chardata(iword)=='snow') gru_struc(iGRU)%hruInfo(iHRU)%nSnow = gru_struc(iGRU)%hruInfo(iHRU)%nSnow+1
               if(chardata(iword)=='soil') gru_struc(iGRU)%hruInfo(iHRU)%nSoil = gru_struc(iGRU)%hruInfo(iHRU)%nSoil+1
               if(chardata(iword)=='snow' .or. chardata(iword)=='soil') exit ! exit once read the layer type
            end do
            deallocate(chardata)
         end do
         deallocate(dataLines)
         
      end do  ! looping through HRUs
   end do  ! looping through HRUs

! **** END OF TEMPORARY CODE ***

! *****************************************************************************
! (3b) allocate space for data structures
! *****************************************************************************

! loop through data structures
   print*, size(structInfo)
   do iStruct=1,size(structInfo)
      print*, iStruct,trim(structInfo(iStruct)%structName) 
      ! allocate space
      select case(trim(structInfo(iStruct)%structName))
      case('time'); call allocGlobal(time_meta,  summa1_struc(n)%timeStruct,  err, message)   ! model forcing data
      case('forc'); call allocGlobal(forc_meta,  summa1_struc(n)%forcStruct,  err, message)   ! model forcing data
      case('attr'); call allocGlobal(attr_meta,  summa1_struc(n)%attrStruct,  err, message)   ! local attributes for each HRU
      case('type'); call allocGlobal(type_meta,  summa1_struc(n)%typeStruct,  err, message)   ! local classification of soil veg etc. for each HRU
      case('mpar'); call allocGlobal(mpar_meta,  summa1_struc(n)%mparStruct,  err, message)   ! model parameters
      case('indx'); call allocGlobal(indx_meta,  summa1_struc(n)%indxStruct,  err, message)   ! model variables
      case('prog'); call allocGlobal(prog_meta,  summa1_struc(n)%progStruct,  err, message)   ! model prognostic (state) variables
      case('diag'); call allocGlobal(diag_meta,  summa1_struc(n)%diagStruct,  err, message)   ! model diagnostic variables
      case('flux'); call allocGlobal(flux_meta,  summa1_struc(n)%fluxStruct,  err, message)   ! model fluxes
      case('bpar'); call allocGlobal(bpar_meta,  summa1_struc(n)%bparStruct,  err, message)   ! basin-average parameters
      case('bvar'); call allocGlobal(bvar_meta,  summa1_struc(n)%bvarStruct,  err, message)   ! basin-average variables
      case('deriv');cycle   ! model derivatives -- no need to have the GRU dimension
      case default; print*, 'unable to identify lookup structure'
      end select

   enddo  ! looping through data structures
   
! allocate space for default model parameters
! NOTE: This is done here, rather than in the loop above, because dpar is not one of the "standard" data structures
   call allocGlobal(mpar_meta,  summa1_struc(n)%dparStruct,  err, message)   ! default model parameters
   allocate(summa1_struc(n)%dt_init(summa1_struc(n)%nGRU),&
        summa1_struc(n)%upArea(summa1_struc(n)%nGRU),&
        summa1_struc(n)%computeVegFlux(summa1_struc(n)%nGRU),stat=err)

! allocate space for the HRUs
   do iGRU=1,summa1_struc(n)%nGRU
      hruCount = gru_struc(iGRU)%hruCount
      allocate(summa1_struc(n)%dt_init(iGRU)%hru(hruCount),&
           summa1_struc(n)%upArea(iGRU)%hru(hruCount),&
           summa1_struc(n)%computeVegFlux(iGRU)%hru(hruCount),stat=err)
!      call LIS_verify(err,'problem allocating space for dt_init, upArea, or computeVegFlux [HRU]')
   end do

! read local attributes for each HRU
   call read_attrb(summa1_struc(n)%nGRU,summa1_struc(n)%nHRU,summa1_struc(n)%attrStruct,&
        summa1_struc(n)%typeStruct,err,message)!; call LIS_verify(err,message)

! *****************************************************************************
! (4a) read description of model forcing datafile used in each HRU
! *****************************************************************************
   call ffile_info(summa1_struc(n)%nHRU,err,message)!; call LIS_verify(err,message)

! *****************************************************************************
! (4b) read model decisions
! *****************************************************************************
   call mDecisions(err,message)!; call LIS_verify(err,message)

! *****************************************************************************
! (3c) allocate space for output statistics data structures
! *****************************************************************************
! loop through data structures
   do iStruct=1,size(structInfo)
 ! allocate space
      select case(trim(structInfo(iStruct)%structName))
      case('forc'); call allocStat(forc_meta, summa1_struc(n)%forcStat,err,message)   ! model forcing data
      case('prog'); call allocStat(prog_meta, summa1_struc(n)%progStat,err,message)   ! model prognostic (state) variables
      case('diag'); call allocStat(diag_meta, summa1_struc(n)%diagStat,err,message)   ! model diagnostic variables
      case('flux'); call allocStat(flux_meta, summa1_struc(n)%fluxStat,err,message)   ! model fluxes
      case('indx'); call allocStat(flux_meta, summa1_struc(n)%indxStat,err,message)   ! index vars
      case('bvar'); call allocStat(bvar_meta, summa1_struc(n)%bvarStat,err,message)   ! basin-average variables
      case default; cycle;
      endselect
 ! check errors
!      call LIS_verify(err,trim(message)//'[statistics for =  '//trim(structInfo(iStruct)%structName)//']')
   enddo ! iStruct

! *****************************************************************************
! (5a) read default model parameters
! *****************************************************************************
   ! read default values and constraints for model parameters (local column, and basin-average)
   call read_pinit(LOCALPARAM_INFO,.TRUE., mpar_meta,localParFallback,err,message)!; call LIS_verify(err,message)
   call read_pinit(BASINPARAM_INFO,.FALSE.,bpar_meta,basinParFallback,err,message)!; call LIS_verify(err,message)
   

! *****************************************************************************
! (5b) read Noah vegetation and soil tables
! *****************************************************************************

! define monthly fraction of green vegetation
!                           J        F        M        A        M        J        J        A        S        O        N        D
   summa1_struc(n)%greenVegFrac_monthly = (/0.01_dp, 0.02_dp, 0.03_dp, 0.07_dp, 0.50_dp, 0.90_dp, 0.95_dp, 0.96_dp, 0.65_dp, 0.24_dp, 0.11_dp, 0.02_dp/)

! read Noah soil and vegetation tables
   call soil_veg_gen_parm(trim(SETNGS_PATH)//'VEGPARM.TBL',                              & ! filename for vegetation table
        trim(SETNGS_PATH)//'SOILPARM.TBL',                             & ! filename for soils table
        trim(SETNGS_PATH)//'GENPARM.TBL',                              & ! filename for general table
        trim(model_decisions(iLookDECISIONS%vegeParTbl)%cDecision),    & ! classification system used for vegetation
        trim(model_decisions(iLookDECISIONS%soilCatTbl)%cDecision))      ! classification system used for soils

! read Noah-MP vegetation tables
   call read_mp_veg_parameters(trim(SETNGS_PATH)//'MPTABLE.TBL',                         & ! filename for Noah-MP table
        trim(model_decisions(iLookDECISIONS%vegeParTbl)%cDecision)) ! classification system used for vegetation
   
   ! define urban vegetation category
   select case(trim(model_decisions(iLookDECISIONS%vegeParTbl)%cDecision))
   case('USGS');                     urbanVegCategory =    1
   case('MODIFIED_IGBP_MODIS_NOAH'); urbanVegCategory =   13
   case('plumberCABLE');             urbanVegCategory = -999
   case('plumberCHTESSEL');          urbanVegCategory = -999
   case('plumberSUMMA');             urbanVegCategory = -999
   case default!; call LIS_verify(30,'unable to identify vegetation category')
   end select

! set default model parameters
   do iGRU=1,summa1_struc(n)%nGRU
      do iHRU=1,gru_struc(iGRU)%hruCount
  ! set parmameters to their default value
         summa1_struc(n)%dparStruct%gru(iGRU)%hru(iHRU)%var(:) = &
              localParFallback(:)%default_val         ! x%hru(:)%var(:)
  ! overwrite default model parameters with information from the Noah-MP tables
         call pOverwrite(&
              summa1_struc(n)%typeStruct%gru(iGRU)%hru(iHRU)%var(iLookTYPE%vegTypeIndex),  &  ! vegetation category
              summa1_struc(n)%typeStruct%gru(iGRU)%hru(iHRU)%var(iLookTYPE%soilTypeIndex), &  ! soil category
              summa1_struc(n)%dparStruct%gru(iGRU)%hru(iHRU)%var,                          &  ! default model parameters
              err,message)!; call LIS_verify(err,message)            ! error control
  ! copy over to the parameter structure
         summa1_struc(n)%mparStruct%gru(iGRU)%hru(iHRU)%var(:) = &
              summa1_struc(n)%dparStruct%gru(iGRU)%hru(iHRU)%var(:)
      end do  ! looping through HRUs
 ! set default for basin-average parameters
      summa1_struc(n)%bparStruct%gru(iGRU)%var(:) = basinParFallback(:)%default_val
   end do  ! looping through GRUs

! *****************************************************************************
! (5c) read trial model parameter values for each HRU, and populate initial data structures
! *****************************************************************************
   call read_param(summa1_struc(n)%nHRU,summa1_struc(n)%typeStruct,&
        summa1_struc(n)%mparStruct,err,message)!; call LIS_verify(err,message)

! *****************************************************************************
! (5d) compute derived model variables that are pretty much constant for the basin as a whole
! *****************************************************************************

! loop through GRUs
   do iGRU=1,summa1_struc(n)%nGRU

 ! calculate the fraction of runoff in future time steps
      call fracFuture(summa1_struc(n)%bparStruct%gru(iGRU)%var,    &  ! vector of basin-average model parameters
           summa1_struc(n)%bvarStruct%gru(iGRU),        &  ! data structure of basin-average variables
           err,message)                    ! error control
! call LIS_verify(err,message)

 ! loop through local HRUs
      do iHRU=1,gru_struc(iGRU)%hruCount
         
         kHRU=0
  ! check the network topology (only expect there to be one downslope HRU)
         do jHRU=1,summa1_struc(n)%nHRU
            if(summa1_struc(n)%typeStruct%gru(iGRU)%hru(iHRU)%var(iLookTYPE%downHRUindex) == &
                 summa1_struc(n)%typeStruct%gru(iGRU)%hru(jHRU)%var(iLookTYPE%hruIndex))then
               if(kHRU==0)then  ! check there is a unique match
                  kHRU=jHRU
               else
!             call LIS_verify(20,'multi_driver: only expect there to be one downslope HRU')
               endif  ! (check there is a unique match)
            endif  ! (if identified a downslope HRU)
         end do
  
  ! check that the parameters are consistent
         call paramCheck(summa1_struc(n)%mparStruct%gru(iGRU)%hru(iHRU)%var,err,message)!; call LIS_verify(err,message)
  
  ! calculate a look-up table for the temperature-enthalpy conversion 
         call E2T_lookup(summa1_struc(n)%mparStruct%gru(iGRU)%hru(iHRU)%var,err,message)!; call LIS_verify(err,message)
 
  ! read description of model initial conditions -- also initializes model structure components
  ! NOTE: at this stage the same initial conditions are used for all HRUs -- need to modify
         call read_icond(gru_struc(iGRU)%hruInfo(iHRU)%nSnow, & ! number of snow layers
              gru_struc(iGRU)%hruInfo(iHRU)%nSoil, & ! number of soil layers
              summa1_struc(n)%mparStruct%gru(iGRU)%hru(iHRU)%var,  & ! vector of model parameters
              summa1_struc(n)%indxStruct%gru(iGRU)%hru(iHRU),      & ! data structure of model indices
              summa1_struc(n)%progStruct%gru(iGRU)%hru(iHRU),      & ! model prognostic (state) variables
              err,message)                           ! error control
         !    call LIS_verify(err,message)
  
  ! re-calculate height of each layer
         call calcHeight(&
              ! input/output: data structures
              summa1_struc(n)%indxStruct%gru(iGRU)%hru(iHRU),   & ! intent(in): layer type
              summa1_struc(n)%progStruct%gru(iGRU)%hru(iHRU),   & ! intent(inout): model prognostic (state) variables for a local HRU
                  ! output: error control
              err,message)!; call LIS_verify(err,message)
  
  ! calculate vertical distribution of root density
         call rootDensty(&
              summa1_struc(n)%mparStruct%gru(iGRU)%hru(iHRU)%var,& ! vector of model parameters
              summa1_struc(n)%indxStruct%gru(iGRU)%hru(iHRU),    & ! data structure of model indices
              summa1_struc(n)%progStruct%gru(iGRU)%hru(iHRU),    & ! data structure of model prognostic (state) variables
              summa1_struc(n)%diagStruct%gru(iGRU)%hru(iHRU),    & ! data structure of model diagnostic variables
              err,message)                         ! error control
!    call LIS_verify(err,message) 
  
  ! calculate saturated hydraulic conductivity in each soil layer
         call satHydCond(summa1_struc(n)%mparStruct%gru(iGRU)%hru(iHRU)%var,& ! vector of model parameters
              summa1_struc(n)%indxStruct%gru(iGRU)%hru(iHRU),    & ! data structure of model indices
              summa1_struc(n)%progStruct%gru(iGRU)%hru(iHRU),    & ! data structure of model prognostic (state) variables
              summa1_struc(n)%fluxStruct%gru(iGRU)%hru(iHRU),    & ! data structure of model fluxes 
              err,message)                         ! error control
!    call LIS_verify(err,message)
         
  ! calculate "short-cut" variables such as volumetric heat capacity
         call v_shortcut(summa1_struc(n)%mparStruct%gru(iGRU)%hru(iHRU)%var,& ! vector of model parameters
              summa1_struc(n)%diagStruct%gru(iGRU)%hru(iHRU),    & ! data structure of model diagnostic variables
              err,message)                         ! error control
!    call LIS_verify(err,message)
  
  ! overwrite the vegetation height
         HVT(summa1_struc(n)%typeStruct%gru(iGRU)%hru(iHRU)%var(iLookTYPE%vegTypeIndex)) = &
              summa1_struc(n)%mparStruct%gru(iGRU)%hru(iHRU)%var(iLookPARAM%heightCanopyTop)
         HVB(summa1_struc(n)%typeStruct%gru(iGRU)%hru(iHRU)%var(iLookTYPE%vegTypeIndex)) = &
              summa1_struc(n)%mparStruct%gru(iGRU)%hru(iHRU)%var(iLookPARAM%heightCanopyBottom)
         
  ! overwrite the tables for LAI and SAI
         if(model_decisions(iLookDECISIONS%LAI_method)%iDecision == specified)then
            SAIM(summa1_struc(n)%typeStruct%gru(iGRU)%hru(iHRU)%var(iLookTYPE%vegTypeIndex),:) = &
                 summa1_struc(n)%mparStruct%gru(iGRU)%hru(iHRU)%var(iLookPARAM%winterSAI)
            LAIM(summa1_struc(n)%typeStruct%gru(iGRU)%hru(iHRU)%var(iLookTYPE%vegTypeIndex),:) = &
                 summa1_struc(n)%mparStruct%gru(iGRU)%hru(iHRU)%var(iLookPARAM%summerLAI)*summa1_struc(n)%greenVegFrac_monthly
         endif
  
  ! initialize canopy drip
  ! NOTE: canopy drip from the previous time step is used to compute throughfall for the current time step
         summa1_struc(n)%fluxStruct%gru(iGRU)%hru(iHRU)%var(iLookFLUX%scalarCanopyLiqDrainage)%dat(1) = 0._dp  ! not used
  
      end do  ! (looping through HRUs)

 ! compute total area of the upstream HRUS that flow into each HRU
      do iHRU=1,gru_struc(iGRU)%hruCount
         summa1_struc(n)%upArea(iGRU)%hru(iHRU) = 0._dp
         do jHRU=1,summa1_struc(n)%nHRU
   ! check if jHRU flows into iHRU
            if(summa1_struc(n)%typeStruct%gru(iGRU)%hru(jHRU)%var(iLookTYPE%downHRUindex)==&
                 summa1_struc(n)%typeStruct%gru(iGRU)%hru(iHRU)%var(iLookTYPE%hruIndex))then
               summa1_struc(n)%upArea(iGRU)%hru(iHRU) = summa1_struc(n)%upArea(iGRU)%hru(iHRU) + &
                    summa1_struc(n)%attrStruct%gru(iGRU)%hru(jHRU)%var(iLookATTR%HRUarea)
            endif   ! (if jHRU is an upstream HRU)
         end do  ! jHRU
      end do  ! iHRU

 ! identify the total basin area for a GRU (m2)
      associate(totalArea => summa1_struc(n)%bvarStruct%gru(iGRU)%var(iLookBVAR%basin__totalArea)%dat(1) )
      totalArea = 0._dp
      do iHRU=1,gru_struc(iGRU)%hruCount
         totalArea = totalArea + summa1_struc(n)%attrStruct%gru(iGRU)%hru(iHRU)%var(iLookATTR%HRUarea)
      end do
      end associate

 ! initialize aquifer storage
 ! NOTE: this is ugly: need to add capabilities to initialize basin-wide state variables
 ! There are two options for groundwater:
 !  (1) where groundwater is included in the local column (i.e., the HRUs); and
 !  (2) where groundwater is included for the single basin (i.e., the GRUS, where multiple HRUS drain into a GRU).
 ! For water balance calculations it is important to ensure that the local aquifer storage is zero if groundwater is treated as a basin-average state variable (singleBasin);
 !  and ensure that basin-average aquifer storage is zero when groundwater is included in the local columns (localColumn).
      select case(model_decisions(iLookDECISIONS%spatial_gw)%iDecision)
  ! the basin-average aquifer storage is not used if the groundwater is included in the local column
      case(localColumn)
         summa1_struc(n)%bvarStruct%gru(iGRU)%var(iLookBVAR%basin__AquiferStorage)%dat(1) = 0._dp ! set to zero to be clear that there is no basin-average aquifer storage in this configuration 
  ! NOTE: the local column aquifer storage is not used if the groundwater is basin-average
  ! (i.e., where multiple HRUs drain to a basin-average aquifer)
      case(singleBasin)
         summa1_struc(n)%bvarStruct%gru(iGRU)%var(iLookBVAR%basin__AquiferStorage)%dat(1) = 1._dp
         do iHRU=1,gru_struc(iGRU)%hruCount
            summa1_struc(n)%progStruct%gru(iGRU)%hru(iHRU)%var(iLookPROG%scalarAquiferStorage)%dat(1) = 0._dp  ! set to zero to be clear that there is no local aquifer storage in this configuration
         end do
      case default!; call LIS_verify(20,'unable to identify decision for regional representation of groundwater')
      end select

 ! initialize time step length for each HRU
      do iHRU=1,summa1_struc(n)%nHRU
         summa1_struc(n)%dt_init(iGRU)%hru(iHRU) = summa1_struc(n)%progStruct%gru(iGRU)%hru(iHRU)%var(iLookPROG%dt_init)%dat(1) ! seconds
      end do
      
   end do  ! (looping through GRUs)

! *****************************************************************************
! (5e) initialize first output sequence 
! *****************************************************************************
! define the file
! NOTE: currently assumes that nSoil is constant across the model domain
   write(fileout,'(a,i0,a,i0,a)') trim(OUTPUT_PATH)//trim(OUTPUT_PREFIX)//'_spinup'//trim(output_fileSuffix)
   call def_output(summa1_struc(n)%nHRU,gru_struc(1)%hruInfo(1)%nSoil,fileout,err,message)!; call LIS_verify(err,message)
   
! write local model attributes and parameters to the model output file
   do iGRU=1,summa1_struc(n)%nGRU
      do iHRU=1,gru_struc(iGRU)%hruCount
         call writeParm(iHRU,summa1_struc(n)%attrStruct%gru(iGRU)%hru(iHRU)%var,attr_meta,err,message)!; call LIS_verify(err,message)
         call writeParm(iHRU,summa1_struc(n)%typeStruct%gru(iGRU)%hru(iHRU)%var,type_meta,err,message)!; call LIS_verify(err,message)
         call writeParm(iHRU,summa1_struc(n)%mparStruct%gru(iGRU)%hru(iHRU)%var,mpar_meta,err,message)!; call LIS_verify(err,message)
      enddo ! HRU
      call writeParm(integerMissing,summa1_struc(n)%bparStruct%gru(iGRU)%var,bpar_meta,err,message)!; call LIS_verify(err,message)
   enddo ! GRU
#endif
enddo
end subroutine summa1_lsm_ini

 ! **************************************************************************************************
 ! private subroutine SOIL_VEG_GEN_PARM: Read soil, vegetation and other model parameters (from NOAH)
 ! **************************************************************************************************
!-----------------------------------------------------------------
SUBROUTINE SOIL_VEG_GEN_PARM(FILENAME_VEGTABLE, FILENAME_SOILTABLE, FILENAME_GENERAL, MMINLU, MMINSL)
!-----------------------------------------------------------------
  use module_sf_noahlsm, only : shdtbl, nrotbl, rstbl, rgltbl, &
       &                        hstbl, snuptbl, maxalb, laimintbl, &
       &                        bb, drysmc, f11, maxsmc, laimaxtbl, &
       &                        emissmintbl, emissmaxtbl, albedomintbl, &
       &                        albedomaxtbl, wltsmc, qtz, refsmc, &
       &                        z0mintbl, z0maxtbl, &
       &                        satpsi, satdk, satdw, &
       &                        theta_res, theta_sat, vGn_alpha, vGn_n, k_soil, &  ! MPC add van Genutchen parameters
       &                        fxexp_data, lvcoef_data, &
       &                        lutype, maxalb, &
       &                        slope_data, frzk_data, bare, cmcmax_data, &
       &                        cfactr_data, csoil_data, czil_data, &
       &                        refkdt_data, natural, refdk_data, &
       &                        rsmax_data, salp_data, sbeta_data, &
       &                        zbot_data, smhigh_data, smlow_data, &
       &                        lucats, topt_data, slcats, slpcats, sltype

  IMPLICIT NONE

  CHARACTER(LEN=*), INTENT(IN) :: FILENAME_VEGTABLE, FILENAME_SOILTABLE, FILENAME_GENERAL
  CHARACTER(LEN=*), INTENT(IN) :: MMINLU, MMINSL
  integer :: LUMATCH, IINDEX, LC, NUM_SLOPE
  integer :: ierr
  INTEGER , PARAMETER :: OPEN_OK = 0

  character*128 :: mess , message

!-----SPECIFY VEGETATION RELATED CHARACTERISTICS :
!             ALBBCK: SFC albedo (in percentage)
!                 Z0: Roughness length (m)
!             SHDFAC: Green vegetation fraction (in percentage)
!  Note: The ALBEDO, Z0, and SHDFAC values read from the following table
!          ALBEDO, amd Z0 are specified in LAND-USE TABLE; and SHDFAC is
!          the monthly green vegetation data
!             CMXTBL: MAX CNPY Capacity (m)
!             NROTBL: Rooting depth (layer)
!              RSMIN: Mimimum stomatal resistance (s m-1)
!              RSMAX: Max. stomatal resistance (s m-1)
!                RGL: Parameters used in radiation stress function
!                 HS: Parameter used in vapor pressure deficit functio
!               TOPT: Optimum transpiration air temperature. (K)
!             CMCMAX: Maximum canopy water capacity
!             CFACTR: Parameter used in the canopy inteception calculati
!               SNUP: Threshold snow depth (in water equivalent m) that
!                     implies 100% snow cover
!                LAI: Leaf area index (dimensionless)
!             MAXALB: Upper bound on maximum albedo over deep snow
!
!-----READ IN VEGETAION PROPERTIES FROM VEGPARM.TBL
!

  OPEN(19, FILE=trim(FILENAME_VEGTABLE),FORM='FORMATTED',STATUS='OLD',IOSTAT=ierr)
  IF(ierr .NE. OPEN_OK ) THEN
     WRITE(message,FMT='(A)') &
          'module_sf_noahlsm.F: soil_veg_gen_parm: failure opening VEGPARM.TBL'
     CALL wrf_error_fatal ( message )
  END IF


  LUMATCH=0

  FIND_LUTYPE : DO WHILE (LUMATCH == 0)
     READ (19,*,END=2002)
     READ (19,*,END=2002)LUTYPE
     READ (19,*)LUCATS,IINDEX

     IF(LUTYPE.EQ.MMINLU)THEN
        WRITE( mess , * ) 'LANDUSE TYPE = ' // TRIM ( LUTYPE ) // ' FOUND', LUCATS,' CATEGORIES'
        ! CALL wrf_message( mess )
        LUMATCH=1
     ELSE
        call wrf_message ( "Skipping over LUTYPE = " // TRIM ( LUTYPE ) )
        DO LC = 1, LUCATS+12
           read(19,*)
        ENDDO
     ENDIF
  ENDDO FIND_LUTYPE
! prevent possible array overwrite, Bill Bovermann, IBM, May 6, 2008
  IF ( SIZE(SHDTBL)       < LUCATS .OR. &
       SIZE(NROTBL)       < LUCATS .OR. &
       SIZE(RSTBL)        < LUCATS .OR. &
       SIZE(RGLTBL)       < LUCATS .OR. &
       SIZE(HSTBL)        < LUCATS .OR. &
       SIZE(SNUPTBL)      < LUCATS .OR. &
       SIZE(MAXALB)       < LUCATS .OR. &
       SIZE(LAIMINTBL)    < LUCATS .OR. &
       SIZE(LAIMAXTBL)    < LUCATS .OR. &
       SIZE(Z0MINTBL)     < LUCATS .OR. &
       SIZE(Z0MAXTBL)     < LUCATS .OR. &
       SIZE(ALBEDOMINTBL) < LUCATS .OR. &
       SIZE(ALBEDOMAXTBL) < LUCATS .OR. &
       SIZE(EMISSMINTBL ) < LUCATS .OR. &
       SIZE(EMISSMAXTBL ) < LUCATS ) THEN
     CALL wrf_error_fatal('Table sizes too small for value of LUCATS in module_sf_noahdrv.F')
  ENDIF

  IF(LUTYPE.EQ.MMINLU)THEN
     DO LC=1,LUCATS
        READ (19,*)IINDEX,SHDTBL(LC),                        &
             NROTBL(LC),RSTBL(LC),RGLTBL(LC),HSTBL(LC), &
             SNUPTBL(LC),MAXALB(LC), LAIMINTBL(LC),     &
             LAIMAXTBL(LC),EMISSMINTBL(LC),             &
             EMISSMAXTBL(LC), ALBEDOMINTBL(LC),         &
             ALBEDOMAXTBL(LC), Z0MINTBL(LC), Z0MAXTBL(LC)
     ENDDO
!
     READ (19,*)
     READ (19,*)TOPT_DATA
     READ (19,*)
     READ (19,*)CMCMAX_DATA
     READ (19,*)
     READ (19,*)CFACTR_DATA
     READ (19,*)
     READ (19,*)RSMAX_DATA
     READ (19,*)
     READ (19,*)BARE
     READ (19,*)
     READ (19,*)NATURAL
  ENDIF
!
2002 CONTINUE

  CLOSE (19)
  IF (LUMATCH == 0) then
     CALL wrf_error_fatal ("Land Use Dataset '"//MMINLU//"' not found in VEGPARM.TBL.")
  ENDIF

!
!-----READ IN SOIL PROPERTIES FROM SOILPARM.TBL
!
  OPEN(19, FILE=trim(FILENAME_SOILTABLE),FORM='FORMATTED',STATUS='OLD',IOSTAT=ierr)
  IF(ierr .NE. OPEN_OK ) THEN
     WRITE(message,FMT='(A)') &
          'module_sf_noahlsm.F: soil_veg_gen_parm: failure opening SOILPARM.TBL'
     CALL wrf_error_fatal ( message )
  END IF

  WRITE(mess,*) 'INPUT SOIL TEXTURE CLASSIFICATION = ', TRIM ( MMINSL )
  ! CALL wrf_message( mess )

  LUMATCH=0



  ! MPC add a new soil table
  FIND_soilTYPE : DO WHILE (LUMATCH == 0)
   READ (19,*)
   READ (19,*,END=2003)SLTYPE
   READ (19,*)SLCATS,IINDEX
   IF(SLTYPE.EQ.MMINSL)THEN
     WRITE( mess , * ) 'SOIL TEXTURE CLASSIFICATION = ', TRIM ( SLTYPE ) , ' FOUND', &
          SLCATS,' CATEGORIES'
     ! CALL wrf_message ( mess )
     LUMATCH=1
   ELSE
    call wrf_message ( "Skipping over SLTYPE = " // TRIM ( SLTYPE ) )
    DO LC = 1, SLCATS
     read(19,*)
    ENDDO
   ENDIF
  ENDDO FIND_soilTYPE
  ! prevent possible array overwrite, Bill Bovermann, IBM, May 6, 2008
  IF ( SIZE(BB    ) < SLCATS .OR. &
       SIZE(DRYSMC) < SLCATS .OR. &
       SIZE(F11   ) < SLCATS .OR. &
       SIZE(MAXSMC) < SLCATS .OR. &
       SIZE(REFSMC) < SLCATS .OR. &
       SIZE(SATPSI) < SLCATS .OR. &
       SIZE(SATDK ) < SLCATS .OR. &
       SIZE(SATDW ) < SLCATS .OR. &
       SIZE(WLTSMC) < SLCATS .OR. &
       SIZE(QTZ   ) < SLCATS  ) THEN
     CALL wrf_error_fatal('Table sizes too small for value of SLCATS in module_sf_noahdrv.F')
  ENDIF

  ! MPC add new soil table
  select case(trim(SLTYPE))
   case('STAS','STAS-RUC')  ! original soil tables
     DO LC=1,SLCATS
        READ (19,*) IINDEX,BB(LC),DRYSMC(LC),F11(LC),MAXSMC(LC),&
             REFSMC(LC),SATPSI(LC),SATDK(LC), SATDW(LC),   &
             WLTSMC(LC), QTZ(LC)
     ENDDO
   case('ROSETTA')          ! new soil table
     DO LC=1,SLCATS
        READ (19,*) IINDEX,&
             ! new soil parameters (from Rosetta)
             theta_res(LC), theta_sat(LC),        &
             vGn_alpha(LC), vGn_n(LC), k_soil(LC), &
             ! original soil parameters
             BB(LC),DRYSMC(LC),F11(LC),MAXSMC(LC),&
             REFSMC(LC),SATPSI(LC),SATDK(LC), SATDW(LC),   &
             WLTSMC(LC), QTZ(LC)
     ENDDO
   case default
     CALL wrf_message( 'SOIL TEXTURE IN INPUT FILE DOES NOT ' )
     CALL wrf_message( 'MATCH SOILPARM TABLE'                 )
     CALL wrf_error_fatal ( 'INCONSISTENT OR MISSING SOILPARM FILE' )
  end select

2003 CONTINUE

  CLOSE (19)

  IF(LUMATCH.EQ.0)THEN
     CALL wrf_message( 'SOIL TEXTURE IN INPUT FILE DOES NOT ' )
     CALL wrf_message( 'MATCH SOILPARM TABLE'                 )
     CALL wrf_error_fatal ( 'INCONSISTENT OR MISSING SOILPARM FILE' )
  ENDIF

!
!-----READ IN GENERAL PARAMETERS FROM GENPARM.TBL
!
  OPEN(19, FILE=trim(FILENAME_GENERAL),FORM='FORMATTED',STATUS='OLD',IOSTAT=ierr)
  IF(ierr .NE. OPEN_OK ) THEN
     WRITE(message,FMT='(A)') &
          'module_sf_noahlsm.F: soil_veg_gen_parm: failure opening GENPARM.TBL'
     CALL wrf_error_fatal ( message )
  END IF

  READ (19,*)
  READ (19,*)
  READ (19,*) NUM_SLOPE

  SLPCATS=NUM_SLOPE
! prevent possible array overwrite, Bill Bovermann, IBM, May 6, 2008
  IF ( SIZE(slope_data) < NUM_SLOPE ) THEN
     CALL wrf_error_fatal('NUM_SLOPE too large for slope_data array in module_sf_noahdrv')
  ENDIF

  DO LC=1,SLPCATS
     READ (19,*)SLOPE_DATA(LC)
  ENDDO

  READ (19,*)
  READ (19,*)SBETA_DATA
  READ (19,*)
  READ (19,*)FXEXP_DATA
  READ (19,*)
  READ (19,*)CSOIL_DATA
  READ (19,*)
  READ (19,*)SALP_DATA
  READ (19,*)
  READ (19,*)REFDK_DATA
  READ (19,*)
  READ (19,*)REFKDT_DATA
  READ (19,*)
  READ (19,*)FRZK_DATA
  READ (19,*)
  READ (19,*)ZBOT_DATA
  READ (19,*)
  READ (19,*)CZIL_DATA
  READ (19,*)
  READ (19,*)SMLOW_DATA
  READ (19,*)
  READ (19,*)SMHIGH_DATA
  READ (19,*)
  READ (19,*)LVCOEF_DATA
  CLOSE (19)

!-----------------------------------------------------------------
END SUBROUTINE SOIL_VEG_GEN_PARM
!-----------------------------------------------------------------
end module summa1_lsmMod

