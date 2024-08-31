!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
! 
! !MODULE: LVT_PRIV_rcMod
! \label(LVT_PRIV_rcMod)
!
! !INTERFACE:
module LVT_PRIV_rcMod 
! 
! !USES:   
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION:
!  
!  Module for specifying model independent variables in LIS 
!  verification toolkit. 
!!
! The variables specified in this module include: 
!
! \begin{description}
!  \item[wout]
!   output format option used in LIS (tiled/gridded/grib/netcdf/grib2) 
!  \item[pass]
!   Number of passes through the time series required for the computation
!   of the metric (e.g. std, anomaly correlations)
!  \item[wtsout]
!   Flag to check whether to output time series data
!  \item[obsCountThreshold]
!   Minimum number of observations (in time), used
!   as a threshold to compute statistics
!  \item[scCountThreshold]
!   Minimum number of points (in time), used
!   as a threshold to compute average seasonal cycles
!  \item[adcCountThreshold]
!   Minimum number of points (in time), used
!   as a threshold to compute average diurnal cycles.
! \end{description}
! !DESCRIPTION: 
! 
! !FILES USED:
!
! !REVISION HISTORY:
!  02 Oct 2008; Sujay Kumar; Initial Specification
!  03 Dec 2012; Shugong Wang; Add lis_version for backward support of LIS 6 
!
!EOP

  implicit none
  type lvtrcdec
     character*500          :: configfile
     integer                :: max_model_types

     character*500           :: runmode
     character*500           :: domain
     integer                 :: nnest
     integer                 :: nDatastreams
     ! The following line is uncommented by Shugong Wang for backword support. 12/04/2012
     character*500           :: lsm
!     character*500           :: routing_model 
!     character*500           :: rtm_model 
     integer                :: nt

     integer                :: metric_sindex
     integer                :: metric_eindex

     integer                :: lsm_index 
     integer                :: lake_index 
     integer                :: glacier_index 
     integer                :: wetland_index 
     integer                :: openwater_index 

     integer                :: npts
     integer                :: ntiles
     integer                :: glbntiles
     integer                :: glbntiles_red
     integer                :: ngrid 
     integer                :: glbngrid
     integer                :: glbngrid_red


     integer                :: gnc
     integer                :: gnr
     integer                :: lnc
     integer                :: lnr
     integer                :: lnc_red
     integer                :: lnr_red
     integer,  allocatable      :: lnc_sc(:)
     integer,  allocatable      :: lnr_sc(:)

     integer                :: pnc  
     integer                :: pnr
     integer                :: input_lnc
     integer                :: input_lnr

     integer                :: nparam

     character*500           :: lvt_out_format
     character*500           :: lvt_wopt
     character*500           :: startmode
     character*500          :: odir
!     character*500           :: model_name
     character*500          :: diagfile

     integer                :: npesx
     integer                :: npesy
     integer                :: halox
     integer                :: haloy

     integer                :: vegsrc

     character*500          :: paramfile
     character*500          :: spTransformAnlysDomain

     integer                :: vfile_form
     real                   :: gridDesc(50)
     real                   :: input_gridDesc(50)

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
     integer                :: doy
     integer                :: yr
     integer                :: mo
     integer                :: prev_mo_rst
     integer                :: prev_mo_tavg
     integer                :: prev_yr_tavg
     integer                :: prev_yr_sout
     integer                :: prev_mo_sout
     integer                :: use_shift_mo
     integer                :: timeAvgOpt
     integer                :: da
     integer                :: hr
     integer                :: mn
     integer                :: ss 
     
     integer                :: nyr
     integer                :: nmo
     integer                :: nda
     integer                :: nhr
     integer                :: nmn
     integer                :: nss

     integer                :: dyr(3)
     integer                :: dmo(3)
     integer                :: dda(3)
     integer                :: dhr(3)
     integer                :: dmn(3)
     integer                :: dss(3)
     integer                :: ddoy(3)
     real*8                 :: dtime(3)      
     real                   :: dgmt(3)

     integer                :: d_nyr(3)
     integer                :: d_nmo(3)
     integer                :: d_nda(3)
     integer                :: d_nhr(3)
     integer                :: d_nmn(3)
     integer                :: d_nss(3)

     integer                :: lyr
     integer                :: lmo
     integer                :: lda
     integer                :: lhr
     integer                :: lmn
     integer                :: lss
     integer                :: ldoy
     real*8                 :: ltime
     real                   :: lgmt

     integer                :: l_nyr
     integer                :: l_nmo
     integer                :: l_nda
     integer                :: l_nhr
     integer                :: l_nmn
     integer                :: l_nss

     real*8                 :: time      
     real                   :: gmt
     
     real*8                 :: etime1
     real                   :: egmt1
     integer                :: edoy1
     integer                :: eyr1
     integer                :: emo1
     integer                :: eda1
     integer                :: ehr1
     integer                :: emn1
     integer                :: ess1
     integer                :: tscount
     integer                :: daycount
     integer                :: daycount_sout
     integer                :: monthCount
     integer                :: monthCount_sout
     character*500          :: rstfile

     integer                :: ts
     character*500           :: tsconv
     integer                :: nts
     integer                :: ntavgs

     integer                :: wtsout
     integer                :: extractts
     integer                :: tavgInterval
     integer, allocatable   :: tlag(:)

     integer                :: wrst
     character*500          :: outputSpecFile
     character*500          :: statsSpecFile
     character*500          :: statsodir
     integer                :: statswriteint
     integer                :: restartInterval
     logical, allocatable   :: resetFlag(:)

     integer                :: nsmlayers
     integer                :: nstlayers
     real                   :: lis_sf_d
     real                   :: lis_rz_d
     integer                :: vinterp_option
     integer :: lis_ts ! EMK

     character*500,  allocatable :: lisvarname(:)
     character*500,  allocatable :: lisvarunit(:)
     character*500,  allocatable :: obsvarname(:,:)
     character*500,  allocatable :: obsvarunit(:,:)

     character*500,  allocatable  :: obssource(:)

     integer                :: smoothObs
     integer                :: pass
     integer                :: curr_pass

     integer                :: dataMask
     integer                :: maskflag
     character*500          :: maskdir
     integer                :: monthly_mask(12)
     integer                :: ftn_summ_file
     integer                :: obsCountThreshold
     integer                :: ComputeErrSC
     integer                :: computeADC
     integer                :: scCountThreshold
     integer                :: adcCountThreshold
     integer                :: scInterval
     integer                :: nasc
     integer                :: nadc
     character*100, allocatable  :: scname(:)
     character*100, allocatable  :: adcname(:)
     logical                :: obs_duplicate
     logical                :: computeFlag
     integer                :: computeEnsMetrics
     integer                :: computeICmetrics
     integer                :: ensLLType
     integer                :: noebal_refet !Reference ET

     real                   :: pval_ci
     integer                :: var_based_strat
     integer                :: var_strat_index
     character*500           :: vname_strat
     real                   :: strat_var_threshold
     integer                :: strat_nlevels

     character*500          :: data_strat_attrib_file
     integer                :: data_based_strat
     character*100, allocatable  :: data_based_strat_file(:)
     character*100, allocatable  :: data_based_strat_var(:)
     integer                :: data_based_nstrats
     real,         allocatable  :: data_based_strat_max(:)
     real,         allocatable  :: data_based_strat_min(:)
     real,         allocatable  :: data_based_strat_delta(:)
     integer,      allocatable  :: data_based_strat_nbins(:)
     real,         allocatable  :: strat_data(:,:)

     integer                :: ntslocs
     character*500          :: tsspecfile
     integer                :: tsspecstyle

     integer                :: n_sc_locs
     character*500          :: sc_specfile
     integer                :: sc_specstyle

     integer                :: n_adc_locs
     character*500          :: adc_specfile
     integer                :: adc_specstyle

     real                   :: udef
     character*500           :: security_class
     character*500           :: distribution_class
     character*500           :: data_category
     character*500           :: area_of_data
     character*500          :: institution = 'NASA GSFC'

     integer                :: nscales
     logical                :: chkTS

     integer                :: anomalyTlength
     character*500          :: sp_avg_mode         
     character*500          :: reg_maskfile
     real,         allocatable  :: regmask(:,:)
     real                   :: regmask_max
     ! The following lines are added by Shugong Wang for backward support of LIS 6 
     integer                :: lis_version
     character*3            :: expcode
     logical                :: lis_output_obs
     
     integer                :: nensem

     real,     allocatable      :: rlat_dn(:)
     real,     allocatable      :: rlon_dn(:)
     real,     allocatable      :: w11_dn(:)
     real,     allocatable      :: w12_dn(:)
     real,     allocatable      :: w21_dn(:)
     real,     allocatable      :: w22_dn(:)
     integer,  allocatable      :: n11_dn(:)
     integer,  allocatable      :: n12_dn(:)
     integer,  allocatable      :: n21_dn(:)
     integer,  allocatable      :: n22_dn(:)

     integer,  allocatable      :: n11_up(:)

     logical                    :: ds1_dup
     logical                    :: ds2_dup
     
     character*500              :: trainingAlg

     integer                :: grib_table
     integer                :: grib_center_id
     integer                :: grib_subcenter_id
     integer                :: grib_grid_id
     integer                :: grib_process_id
     character*50           :: grib_packing_type

     integer                    :: HYCOM_nc
     integer                    :: HYCOM_nr
     ! EMK...Add support for ARC and ANT sea ice fields
     integer                    :: HYCOM_aice_arc_nc
     integer                    :: HYCOM_aice_arc_nr
     integer                    :: HYCOM_aice_ant_nc
     integer                    :: HYCOM_aice_ant_nr
     integer                    :: HYCOM_hi_arc_nc
     integer                    :: HYCOM_hi_arc_nr
     integer                    :: HYCOM_hi_ant_nc
     integer                    :: HYCOM_hi_ant_nr
     logical                    :: HYCOM_proc_start
     real,  allocatable         :: HYCOM_n11(:)
     real,  allocatable         :: HYCOM_aice_arc_n11(:)
     real,  allocatable         :: HYCOM_aice_ant_n11(:)
     real,  allocatable         :: HYCOM_hi_arc_n11(:)
     real,  allocatable         :: HYCOM_hi_ant_n11(:)
     integer                    :: processHYCOM
     character*100              :: HYCOMdir
     integer                    :: applyNoiseReductionFilter
     character*100              :: smoothingFilterType

     ! For USAFSIpost
     character(len=10) :: yyyymmddhh
     logical :: output_native
     logical :: output_global_ll0p25
     logical :: output_nh_ps16
     logical :: output_sh_ps16
     logical :: output_nh_ps16_snodep
     logical :: output_sh_ps16_snodep
     character(len=255) :: input_dir
     character(len=255) :: input_prefix
     character(len=255) :: output_dir
  end type lvtrcdec
  
  type lisrcdec
     integer                :: ts     
     character*500           :: anlys_data_class
     integer                :: nsf_model_types
     integer, allocatable       :: sf_model_type(:)
     character*500, allocatable  :: sf_model_type_name(:)
     integer, allocatable       :: sf_model_type_select(:)
     character*500, allocatable  :: sf_model_type_name_select(:)

     integer                :: nsurfacetypes
     character*500          :: domfile
     character*500          :: map_proj
     character*500          :: odir
     integer                :: nest
     character*500          :: style
     character*500          :: format

     integer                :: bareclass 
     integer                :: urbanclass
     integer                :: snowclass 
     integer                :: waterclass
     integer                :: wetlandclass
     integer                :: glacierclass
     integer                :: nvegtypes
     integer                :: nelevbands
     integer                :: nslopebands
     integer                :: naspectbands
     integer                :: nsoiltypes
     integer                :: nsoilfbands

     integer                :: ntiles
     integer                :: glbntiles
     integer                :: glbntiles_red
     integer                :: ngrid
     integer                :: glbngrid
     integer                :: glbngrid_red
     integer, allocatable   :: npatch(:)
     integer, allocatable   :: glbnpatch(:)
     integer, allocatable   :: glbnpatch_red(:)

     integer                :: lnc
     integer                :: lnr
     integer                :: lnc_red
     integer                :: lnr_red
     integer                :: gnc
     integer                :: gnr
     integer                :: nensem
     character*500           :: wopt
     character*500          :: useelevationmap
     character*500          :: useslopemap
     character*500          :: useaspectmap
     character*500          :: usetexturemap
     character*500          :: usesoilfractionmap
     real                   :: gridDesc(50)

     integer                :: surface_maxt
     real                   :: surface_minp    
     integer                :: soilt_maxt
     real                   :: soilt_minp    
     integer                :: soilf_maxt
     real                   :: soilf_minp    
     integer                :: elev_maxt
     real                   :: elev_minp    
     integer                :: slope_maxt
     real                   :: slope_minp    
     integer                :: aspect_maxt
     real                   :: aspect_minp    
     character*500          :: outputSpecFile

     real,     allocatable      :: rlat_dn(:)
     real,     allocatable      :: rlon_dn(:)
     real,     allocatable      :: w11_dn(:)
     real,     allocatable      :: w12_dn(:)
     real,     allocatable      :: w21_dn(:)
     real,     allocatable      :: w22_dn(:)
     integer,  allocatable      :: n11_dn(:)
     integer,  allocatable      :: n12_dn(:)
     integer,  allocatable      :: n21_dn(:)
     integer,  allocatable      :: n22_dn(:)

     integer,  allocatable      :: n11_up(:)
     character*500               :: model_name
     integer                :: nsmlayers
     integer                :: nstlayers
     real,          allocatable :: smthick(:)
     real,          allocatable :: stthick(:)
     real,          allocatable :: smdepth(:)
     real,          allocatable :: stdepth(:)
     
     real,          allocatable :: vic_depth(:,:)


  end type lisrcdec
end module LVT_PRIV_rcMod
