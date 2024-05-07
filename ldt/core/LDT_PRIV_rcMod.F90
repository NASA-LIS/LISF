!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
module LDT_PRIV_rcMod 
!BOP
!
! !MODULE: LDT_PRIV_rcMod
!
! !DESCRIPTION:
!  
!  Module for specifying model independent variables in Land
!  Data Toolkit (LDT)
!
! !REVISION HISTORY:
!  07 Jul 2010; Sujay Kumar; Initial Specification
!
! The variables specified in this module include: 
!
! \begin{description}
!  \item[wout]
!    output format option used in LDT (tiled/gridded/grib/netcdf/grib2) 
!  \item[pass]
!    Number of passes through the time series required for the computation
!    of the metric (e.g. std, anomaly correlations)
! \end{description}
!
!EOP
  use LDT_constantsMod, only : LDT_CONST_PATH_LEN

  implicit none

  type ldtrcdec

     character*50           :: runmode
     integer                :: nnest
     character(len=LDT_CONST_PATH_LEN), allocatable :: paramAttribsFile(:)
     character*50,  allocatable :: lis_map_proj(:)
     real, allocatable          :: lis_map_resfactor(:)
     integer                :: add_buffer   ! KA
     integer                :: x_buffer
     integer                :: y_buffer
     integer                :: allmaskland  ! KA: Temporary for NoahMP+SM coupling

! -- Land surface input parameters:
     integer                :: max_model_types 
     integer                :: nsf_model_types
     integer, allocatable       :: sf_model_type(:)
     character*50, allocatable  :: sf_model_type_name(:)
     integer, allocatable       :: sf_model_type_select(:)
     character*50, allocatable  :: sf_model_type_name_select(:)

     integer                :: lsm_index 
     integer                :: lake_index 
     integer                :: glacier_index 
     integer                :: wetland_index 
     integer                :: openwater_index 

! -- Land/Surface model inputs:
     character*50           :: lsm
     character*50           :: lakemodel
     character*50           :: routingmodel
     logical                :: inc_water_pts
     real,          allocatable :: gridcell_water_frac(:)
     real,          allocatable :: gridcell_glacier_frac(:)
   ! LSM/Crop-specific entries:
     logical,       allocatable :: assimcropinfo(:)

! -- Forcing input parameters:
     integer                :: nf
     integer                :: nmetforc
     integer                :: nperforc
     integer                :: nmetforc_parms
     character*50           :: metforc_blend_alg
     character*50, allocatable :: metforc(:)
     character*50, allocatable :: met_gridtransform(:)
     character*50, allocatable :: met_ecor(:)
     integer,      allocatable :: met_nf(:)
     real,         allocatable :: met_ts(:)
     character*50, allocatable :: met_tinterp(:)
     logical,      allocatable :: met_zterp(:)
     integer,      allocatable :: met_validhr(:)

     character*50, allocatable :: metforc_parms(:)
     character*50, allocatable :: metforc_parmsrc(:)
     character*50, allocatable :: met_gridtransform_parms(:)
     character*50, allocatable :: met_ecor_parms(:)
     character*50, allocatable :: met_proj(:)
     real,         allocatable :: met_gridDesc(:,:)
     integer,      allocatable :: met_nc(:)
     integer,      allocatable :: met_nr(:)

     integer, allocatable      :: met_nensem(:)
     integer,      allocatable :: pcp_downscale(:)

     real                   :: metForcOutInterval
     real                   :: metForcProcInterval
     integer                :: metForcTWstarthour  ! KRA
     integer                :: metForcTWendhour    ! KRA

     character(len=LDT_CONST_PATH_LEN)          :: forcvarlistFile
     character(len=LDT_CONST_PATH_LEN)          :: forcattribFile
     character(len=LDT_CONST_PATH_LEN)          :: forcpertattribFile

     logical                :: zterp_correction

     integer, allocatable   :: nsf(:)
     integer, allocatable   :: findtime1(:),findtime2(:)

! -- Tile parameters:
     integer                :: nt             ! Number veg types (to be removed)
     integer, allocatable   :: numcrop(:)     ! Number crop types 

     integer, allocatable   :: surface_maxt(:)
     real, allocatable      :: surface_minp(:) 
     integer, allocatable   :: soilt_maxt(:)
     real, allocatable      :: soilt_minp(:)    
     integer, allocatable   :: soilf_maxt(:)
     real, allocatable      :: soilf_minp(:)    
     integer, allocatable   :: elev_maxt(:)
     real, allocatable      :: elev_minp(:)    
     integer, allocatable   :: slope_maxt(:)
     real, allocatable      :: slope_minp(:)    
     integer, allocatable   :: aspect_maxt(:)
     real, allocatable      :: aspect_minp(:)    

     integer, allocatable   :: ntiles(:)
     integer, allocatable   :: glbntiles(:)
     integer, allocatable   :: glbntiles_red(:)
     integer, allocatable   :: npatch(:,:)
     integer, allocatable   :: glbnpatch(:,:)
     integer, allocatable   :: glbnpatch_red(:,:)
     integer, allocatable   :: ngrid(:) 
     integer, allocatable   :: glbngrid(:)
     integer, allocatable   :: glbngrid_red(:)
     integer, allocatable   :: gnc(:)
     integer, allocatable   :: gnr(:)
     integer, allocatable   :: lnc(:)
     integer, allocatable   :: lnr(:)
     integer, allocatable   :: gnc_buf(:)
     integer, allocatable   :: gnr_buf(:)
     integer, allocatable   :: lnc_b(:)
     integer, allocatable   :: lnr_b(:)
     integer, allocatable   :: lnc_red(:)
     integer, allocatable   :: lnr_red(:)
     integer, allocatable   :: ncatg(:)
     integer, allocatable   :: nensem(:)
     integer, allocatable   :: nmaskpts(:)

     logical, allocatable   :: cliplandmask(:)

! -- Time parameters:
     integer                :: doy
     integer                :: yr
     integer                :: mo
     integer                :: da
     integer                :: hr
     integer                :: mn
     integer                :: ss
     integer                :: ms
     real*8                 :: time
     real                   :: gmt
     real                   :: ts
     integer, allocatable   :: tscount(:)
     real, allocatable      :: nts(:)

     integer                :: sdoy        
     integer                :: sss         
     integer                :: smn         
     integer                :: shr         
     integer                :: sda         
     integer                :: smo         
     integer                :: syr         
     integer                :: endcode        
     integer                :: ess            
     integer                :: emn            
     integer                :: edoy           
     integer                :: ehr            
     integer                :: eda            
     integer                :: emo            
     integer                :: eyr            
     integer                :: endtime        
     real*8                 :: etime               
     real                   :: egmt
     real*8                 :: etime1
     real                   :: egmt1
     integer                :: edoy1
     integer                :: eyr1
     integer                :: emo1
     integer                :: eda1
     integer                :: ehr1
     integer                :: emn1
     integer                :: ess1

     integer                :: monthCount
     integer                :: alrm_prev_mo
     
     real                   :: udef
     character*20           :: security_class
     character*20           :: distribution_class
     character*20           :: data_category
     character*20           :: area_of_data

     integer                :: npesx
     integer                :: npesy
     integer                :: halox
     integer                :: haloy

! -- Parameter grid description arrays:
     real, allocatable      :: gridDesc(:,:)

     real, allocatable      :: lc_gridDesc(:,:)
     real, allocatable      :: mask_gridDesc(:,:)
     real, allocatable      :: reg_gridDesc(:,:)
     real, allocatable      :: sfctype_gridDesc(:,:)
     real, allocatable      :: soil_gridDesc(:,:)
     real, allocatable      :: soiltext_gridDesc(:,:)
     real, allocatable      :: hsg_gridDesc(:,:)
     real, allocatable      :: topo_gridDesc(:,:)

     real, allocatable      :: global_mask(:,:)

! -- Parameter projection options:
     character*50           :: lc_proj
     character*50           :: mask_proj
     character*50           :: sfctype_proj
     character*50           :: reg_proj
     character*50           :: soils_proj
     character*50           :: soiltext_proj
     character*50           :: hsg_proj
     character*50           :: topo_proj

! -- Parameter spatial transformation options:
     character*50, allocatable  :: lc_gridtransform(:)
     character*50, allocatable  :: mask_gridtransform(:)
     character*50, allocatable  :: reg_gridtransform(:)
     character*50, allocatable  :: soils_gridtransform(:)
     character*50, allocatable  :: soiltext_gridtransform(:)
     character*50, allocatable  :: topo_gridtransform(:)

     character*50,  allocatable :: lc_type(:)
     character*50,  allocatable :: mask_type(:)
     character*50,  allocatable :: mask_source(:)
     character*50,  allocatable :: soil_classification(:)

     character*50               :: create_soilparms_option

! -- Parameter filepath names:
     character(len=LDT_CONST_PATH_LEN), allocatable :: mfile(:)
     character(len=LDT_CONST_PATH_LEN), allocatable :: vfile(:)
     character(len=LDT_CONST_PATH_LEN), allocatable :: sfctypefile(:)
     character(len=LDT_CONST_PATH_LEN), allocatable :: regfile(:)

     character(len=LDT_CONST_PATH_LEN), allocatable :: safile(:) 
     character(len=LDT_CONST_PATH_LEN), allocatable :: clfile(:) 
     character(len=LDT_CONST_PATH_LEN), allocatable :: sifile(:) 
     character(len=LDT_CONST_PATH_LEN), allocatable :: gravelfile(:) 
     character(len=LDT_CONST_PATH_LEN), allocatable :: txtfile(:)
     character(len=LDT_CONST_PATH_LEN), allocatable :: pofile(:)
     character(len=LDT_CONST_PATH_LEN), allocatable :: psisatfile(:) 
     character(len=LDT_CONST_PATH_LEN), allocatable :: ksatfile(:) 
     character(len=LDT_CONST_PATH_LEN), allocatable :: bexpfile(:) 
     character(len=LDT_CONST_PATH_LEN), allocatable :: qzfile(:)   
     character(len=LDT_CONST_PATH_LEN), allocatable :: dsoilfile(:)
     character(len=LDT_CONST_PATH_LEN), allocatable :: bdrckdepfile(:)
     character(len=LDT_CONST_PATH_LEN), allocatable :: hsgfile(:)
     character(len=LDT_CONST_PATH_LEN), allocatable :: bulkdensfile(:)
     character(len=LDT_CONST_PATH_LEN), allocatable :: domrockfile(:)
     character(len=LDT_CONST_PATH_LEN), allocatable :: rockvolfile(:)
     character(len=LDT_CONST_PATH_LEN), allocatable :: awcfile(:)
     character(len=LDT_CONST_PATH_LEN), allocatable :: permabfile(:)

     character(len=LDT_CONST_PATH_LEN), allocatable :: iscfile(:) 
     character(len=LDT_CONST_PATH_LEN), allocatable :: elevfile(:)
     character(len=LDT_CONST_PATH_LEN), allocatable :: slfile(:)
     character(len=LDT_CONST_PATH_LEN), allocatable :: aspfile(:)
     character(len=LDT_CONST_PATH_LEN), allocatable :: curvfile(:)

     character(len=LDT_CONST_PATH_LEN), allocatable :: glaciermask(:)

   ! GDAS terrain height files:
     character(len=LDT_CONST_PATH_LEN), allocatable :: gdasT126elevfile(:)
     character(len=LDT_CONST_PATH_LEN), allocatable :: gdasT170elevfile(:)
     character(len=LDT_CONST_PATH_LEN), allocatable :: gdasT254elevfile(:)
     character(len=LDT_CONST_PATH_LEN), allocatable :: gdasT382elevfile(:)
     character(len=LDT_CONST_PATH_LEN), allocatable :: gdasT574elevfile(:)
     character(len=LDT_CONST_PATH_LEN), allocatable :: gdasT1534elevfile(:)

   ! ECMWF terrain height files:
     character(len=LDT_CONST_PATH_LEN), allocatable :: elevfileifs23r4(:)  
     character(len=LDT_CONST_PATH_LEN), allocatable :: elevfileifs25r1(:)  
     character(len=LDT_CONST_PATH_LEN), allocatable :: elevfileifs30r1(:)  
     character(len=LDT_CONST_PATH_LEN), allocatable :: elevfileifs33r1(:)  
     character(len=LDT_CONST_PATH_LEN), allocatable :: elevfileifs35r2(:) 
     character(len=LDT_CONST_PATH_LEN), allocatable :: elevfileifs35r3(:)  
     character(len=LDT_CONST_PATH_LEN), allocatable :: elevfileifs36r1(:)  
     character(len=LDT_CONST_PATH_LEN), allocatable :: elevfileifs37r2(:)  

     logical,       allocatable :: monthlyData(:)
     logical,       allocatable :: quarterlyData(:)

! -- Output file attributes:
     integer                :: nobs
     character*50           :: wopt
     character*50           :: wout
     integer                :: wout_form
     integer                :: lis_wopt
     character*3            :: expcode
     character(len=LDT_CONST_PATH_LEN)          :: odir
     character*20           :: model_name
     character(len=LDT_CONST_PATH_LEN)           :: diagfile
     character(len=LDT_CONST_PATH_LEN)          :: mpfillfile
     character(len=LDT_CONST_PATH_LEN), allocatable :: outputSpecFile(:)

     integer                :: wsingle
     character*50           :: wstyle
     integer                :: grib_table
     integer                :: grib_center_id
     integer                :: grib_subcenter_id
     integer                :: grib_grid_id
     integer                :: grib_process_id
     character*50           :: startcode
     integer                :: plevel
     real                   :: tavgInterval
     logical                :: computeFlag

     integer                :: waterclass
     integer                :: lakeclass
     integer                :: bareclass 
     integer                :: urbanclass
     integer                :: snowclass 
     integer                :: wetlandclass
     integer                :: glacierclass
     integer                :: permafrostclass
     integer                :: cropclass1
     integer                :: cropclass2
     integer                :: cropclass3
     integer                :: cropclass4
     integer                :: cropclass5
     integer                :: grassclass
     integer                :: shrubclass1
     integer                :: shrubclass2
  
     integer                :: prev_mo
     character(len=LDT_CONST_PATH_LEN)           :: maskdir

! -- DA Preprocessing Inputs:
     character*50           :: obs_src
     integer                :: comp_cdf
     integer                :: comp_obsGrid
     integer                :: cdf_nbins
     integer                :: cdf_ntimes
     integer                :: group_cdfs
     character*50           :: group_cdfs_attrib_file
     character*50           :: group_cdfs_strat_file
     real                   :: group_cdfs_min
     real                   :: group_cdfs_max
     integer                :: group_cdfs_nbins
     integer                :: daily_interp_switch  !0:on; 1:off (Y.Kwon)

     integer                :: strat_cdfs
     integer                :: write_strat_cdfs
     character*50           :: strat_src
     character*50           :: strat_file
     !character*50           :: strat_cdfs_attrib_file
     integer                :: strat_cdfs_nbins
     integer                :: stratified_cdfs_nbins
     real                   :: strat_cdfs_min
     real                   :: strat_cdfs_max
     integer, allocatable   :: stratification_data(:,:)


     integer                :: sp_sampl_cdfs
     integer                :: sp_sample_cdf_rad

     integer, allocatable   :: cdf_strat_data(:)
     integer                :: anomalyObsProc

     integer                :: obsCountThreshold
     character(len=LDT_CONST_PATH_LEN)          :: dapreprocfile
     integer                :: applyMask
     integer, allocatable   :: obssource(:)
     integer                :: pass
     integer                :: pass_id

     integer                :: ftn_cdf
     integer                :: ftn_strat_cdf
     integer                :: ftn_DAobs_domain
     character*100          :: institution = 'NASA GSFC'     

! -- Metforcing preprocessing Config inputs:
     character*50           :: timeDscaleType

! -- statistical downscaling
     character*50           :: statDscaleType
     character*50           :: statDscaleMode  !downscale/forecast

! -- Restart file parameters:
     integer, allocatable   :: rstflag(:)
     integer, allocatable   :: gridchange(:)
     character*100          :: rstsource
     character*50           :: ensrstmode
     character*50           :: ensrstsampling
     character*140          :: inputrst
     character*140          :: outputrst
     integer                :: nens_in
     integer                :: nens_out

     !ag (1Nov2017)
     integer                :: routing_grid_count

  end type ldtrcdec
  
end module LDT_PRIV_rcMod
