!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
module AGRMET_forcingMod
!BOP
! !MODULE: AGRMET_forcingMod
! 
! !DESCRIPTION: 
!  This module contains the variables and data structures used for 
!  the implementation of Air Force Weather Agency (AFWA)'s Agricultural
!  Meteorological Model (AGRMET) analyses to produce meteorological 
!  forcing variables for land surface modeling. The algorithms use
!  AFWA's unique database of observational and satellite-based 
!  meteorological fields. The input data includes data from Surface
!  Observations database, CDFS II model, and the Special Sensor 
!  Microwave/Imager Sounder (SSMIS) database. External inputs come from 
!  either the NCEP Global Forecast System (GFS), previously known
!  as the Aviation (AVN) model.  AGRMET algorithms also include
!  sophisticated precipitation analysis to compute global precipitation
!  at every grid point around the globe. Estimates of precipitation
!  based on satellite data and climatological tables are blended with
!  surface observations to form a global dataset. 
! 
!  Since the algorithms use a large number of datasets for each analyses,
!  the code assumes the following hierarchy for each instance. For 
!  example, if the root directory for storing the files is 'FORCING/AFWA/'
!  and the current instance is decemeber 1st, 2005, the files are
!  stored under 'FORCING/AFWA/20051201/' directory. Under this directory
!  the files should be stored under the following subdirectories, which
!  can be specified through the lis configuration file. 
!
!  \begin{verbatim}
!   <JMOBS data directory>       : surface and precip 
!                                  obs (sfcobs_* and preobs_*)
!   <GEOPRECIP data directory>   : GEOPRECIP files (prec08* and rank08*) 
!   <surface fields directory>   : processed surface fields (sfc*)
!   <GFS data directory>         : GFS data (MT.avn*)
!   <merged precip directory>    : merged precip
!                                  output (premrgp*)
!   <SSMIS data directory>       : SSMIS data (ssmira_*)
!   <cloud data directory>       : WWMCA data (WWMCA*)
!  \end{verbatim}
!
!  The implementation of AGRMET algorithms use a common derived
!  data type {\tt agrmet\_struc} that includes variables that specify
!  runtime options, information needed for geospatial transformation
!  and other variables. These variables are described below. 
! 
!  \begin{description}
!  \item[radProcessInterval]
!   interval in hours to perform a radiation analysis
!  \item[pcpProcessInterval]
!   interval in hours to perform a precipitation analysis
!  \item[readgfsInterval]
!   interval in hours to read the GFS data
!  \item[radProcessAlarmTime]
!   alarm time to keep track of the radiation processing frequency
!  \item[pcpProcessAlarmTime]
!   alarm time to keep track of the precip processing frequency
!  \item[gfsAlarmTime]
!   alarm time to keep track of the GFS data reading frequency
!  \item[pcpclimoAlarmTime]
!   alarm time to keep track of precip climatology reading frequency
!  \item[albAlarmTime]
!   alarm time to keep track of albedo climatology reading frequency
!  \item[gfracAlarmTime]
!   alarm time to keep track of green fraction climatology reading frequency
!  \item[agrmetdir]
!   root directory for the input files used in different analyses.
!  \item[climodir]
!   precip climatology directory
!  \item[maskfile]
!   landmask file used in the AGRMET analyses 
!  \item[retroFileRoot]
!   Root of filename for input LIS files in retrospective mode
!  \item[terrainfile]
!   terrain file used in the AGRMET analyses 
!  \item[sfcntmfile]
!   file with the spreading radii used for the barnes analysis 
!   on the GFS and surface obs. 
!  \item[analysisdir]
!   directory to write the AGRMET intermediate analyses.
!  \item[agrmettime1]
!    The nearest, previous hourly instance of the incoming 
!    data (as a real time). 
!  \item[agrmettime2]
!    The nearest, next hourly instance of the incoming 
!    data (as a real time).
!  \item[agrmetpcptime1]
!    The nearest, previous hourly instance of the precip 
!    processing (as a real time). 
!  \item[agrmetpcptime2]
!    The nearest, next hourly instance of the precip 
!    processing (as a real time)
!  \item[cmortime]
!    The nearest, hourly instance of the incoming 
!    data (as a real time).
!  \item[pcpobswch]
!    choice for reading precip observations
!  \item[pwswch]
!    choice for using present/past weather estimate
!  \item[raswch]
!    choice for reading SSMIS data
!  \item[cdfs2swch]
!    choice for using CDFS2 estimate
!  \item[clswch]
!    choice for using precip climatology
!  \item[geoswch]
!    choice for using GEOPRECIP data
!  \item[cmorswch]
!    choice for using CMORPH data
!  \item[cmorthere]
!    indicates if cmorph file is available
!  \item[razero]
!    choice for using SSMIS zeros
!  \item[salp]
!    snow distribution shape parameter
!  \item[clmult]
!    alternate monthly weighting factor
!  \item[mnpr]
!    minimum 3 hour climo precip value required to generate
!    a non-zero CDFS2 total cloud based precip estimate
!  \item[mxpr]
!    maxmimum 3 hour climo precip value required to generate
!    a non-zero CDFS2 total cloud based precip estimate
!  \item[mnpd]
!    minimum precip-per-precip day multiplier used to generate
!    a non-zero CDFS2 total cloud based precip estimate
!  \item[mxpd]
!    maxmimum precip-per-precip day multiplier used to generate
!    a non-zero CDFS2 total cloud based precip estimate
!  \item[cldth]
!    cloud threshold to generate a CDFS2 based precipitation
!    estimate
!  \item[mdmv]
!    median cloud cover percentage to move to for the CDFS2
!    based precipitation estimate
!  \item[mdpe]
!    percent to move to median cloud cover percentage 
!    to move to for the CDFS2 based precipitation
!  \item[ovpe]
!    overcast percentage to move to for CDFS2 based 
!    precipitation estimate
!  \item[cdfs2int]
!    CDFS2 time interval to look for cloud amount
!  \item[pcap]
!    3 hour maximum precip ceiling
!  \item[wndwgt]
!    weighting factor for first guess winds
!  \item[minwnd]
!    minimum allowable wind speed on the agrmet grid
!  \item[gridspan]
!    value indicating the span of the LIS domain 
!    (1-NH only, 2-SH only, 3-span includes both NH and SH)
!  \item[mo1]
!    number of elements in the NH for the LIS grid
!  \item[mo2]
!    number of elements in the SH for the LIS grid
!  \item[hemi\_nc]
!    number of columns (in the east-west dimension) for each
!    hemisphere
!  \item[hemi\_nr]
!    number of rows (in the north-south dimension) for each
!    hemisphere   
!  \item[shemi]
!    value indicating the starting index of hemisphere loops
!    (1-NH only,2-SH only, 1-span includes both NH and SH)
!  \item[nhemi]
!    value indicating the ending index of hemisphere loops
!    (1-NH only,2-SH only, 2-span includes both NH and SH)
!  \item[mi]
!    number of elements in the input grid
!  \item[micmor]
!    Number of points in the CMORPH input grid
!  \item[imax]
!    x-dimension size of the input grid
!  \item[jmax]
!    y-dimension size of the input grid
!  \item[interp]
!    variable specifying if spatial interpolation needs to be
!    done
!  \item[land]
!    landmask data
!  \item[veg]
!    vegtype data
!  \item[alt]
!    terrain data
!  \item[irad]
!    radius of observation points for surface OBS
!  \item[radius]
!    radius of observations points for precip OBS
!  \item[rlat1\_nh]
!    latitudes of the input grid in the northern hemisphere
!    to be used for bilinear interpolation
!  \item[rlat1\_sh]
!    latitudes of the input grid in the southern hemisphere
!    to be used for bilinear interpolation
!  \item[rlon1\_nh]
!    longitudes of the input grid in the northern hemisphere
!    to be used for bilinear interpolation
!  \item[rlon1\_sh]
!    longitudes of the input grid in the southern hemisphere
!    to be used for bilinear interpolation
!  \item[n111\_nh,n121\_nh,n211\_nh,n221\_nh]
!    neighbor information of the input grid in the northern 
!    hemisphere to be used for bilinear interpolation
!  \item[n111\_sh,n121\_sh,n211\_sh,n221\_sh]
!    neighbor information of the input grid in the southern
!    hemisphere to be used for bilinear interpolation
!  \item[w111\_nh,w121\_nh,w211\_nh,w221\_nh]
!    bilinear interpolation weights in the northern 
!    hemisphere
!  \item[w111\_sh,w121\_sh,w211\_sh,w221\_sh]
!    bilinear interpolation weights in the southern 
!    hemisphere
!  \item[w112cmor,w122cmor,w212cmor,w222cmor]
!    Arrays containing the weights of the CMORPH grid 
!    for each grid point in LIS, for conservative interpolation.
!  \item[n112\_nh, n112\_sh]
!    neighbor information of the input grid for neighbor
!    interpolation
!  \item[n112cmor,n122cmor,n212cmor,n222cmor]
!    Arrays containing the neighbor information of the CMORPH grid 
!    for each grid point in LIS, for conservative interpolation. 
!  \item[smask1]
!    landmask to used to fill any gaps in bilinear interpolation
!    due to mismatches in the LIS and AGRMET masks
!  \item[smask2]
!    landmask to used to fill any gaps in neighbor interpolation
!    due to mismatches in the LIS and AGRMET masks
!  \item[fillflag1]
!    flag to check for filling gaps due to LIS and AGRMET
!    mask mismatches, for bilinear interpolation
!  \item[fillflag2]
!    flag to check for filling gaps due to LIS and AGRMET
!    mask mismatches, for neighbor interpolation 
!  \item[cliprc1,cliprc2]
!    month 1 and 2 3-hour precip climo values 
!  \item[clippd1,clippd2]
!    month 1 and 2 precip-per-precip day amount
!  \item[clirtn1,clirtn2]
!    month 1 and 2 climatological rtneph percent cloud cover
!  \item[cliprc]
!    3 hour climo precip used for estimate
!  \item[clippd]
!    climatological precip-per-precip day amount
!  \item[clirtn]
!    climatological rtneph percent cloud cover
!  \item[mrgp]
!    merged precip amounts
!  \item[pcp\_start]
!    flag indicating the start of the precip analysis
!  \item[agr\_hgt\_c]
!    geopotential heights on the AGRMET grid for the current time
!  \item[agr\_rh\_c]
!    relative humidity on the AGRMET grid for the current time
!  \item[agr\_tmp\_c]
!    temperature on the AGRMET grid for the current time
!  \item[agr\_wspd\_c]
!    wind speed on the AGRMET grid for the current time
!  \item[agr\_hgt\_p]
!    geopotential heights on the AGRMET grid from the previous 6 hr
!  \item[agr\_rh\_p]
!    relative humidity on the AGRMET grid from the previous 6 hr
!  \item[agr\_tmp\_p]
!    temperature on the AGRMET grid from the previous 6 hr
!  \item[agr\_wspd\_p]
!    wind speed on the AGRMET grid from the previous 6 hr
!  \item[step]
!    time loop counter used to calculate time interpolation weights
!  \item[sfcprs]
!    first guess pressure fields on the AGRMET grid
!  \item[srcrlh]
!    first guess relative humidity fields on the AGRMET grid
!  \item[sfctmp]
!    first guess temperature fields on the AGRMET grid
!  \item[lasprs]
!    first guess pressure fields for ending 6 hour time on the AGRMET grid
!  \item[lasrlh]
!   first guess relative humidity fields for ending 6 hour time
!   on the AGRMET grid
!  \item[lasspd]
!   first guess wind speed fields for ending 6 hour time
!   on the AGRMET grid
!  \item[lastmp]
!   first guess temperature fields for ending 6 hour time
!   on the AGRMET grid
!  \item[max\_sfcobs]
!   size of arrays holding surface observation data
!  \item[max\_pcpobs]
!   size of arrays holding precipitation observation data
!  \item[findtime1, findtime2]
!   boolean flags to indicate which time is to be read for 
!   temporal interpolation.
!  \item[metdata1, metdata2]
!   contains the previous and next forcing values, respectively
!  \end{description}
!
! !USES:
  use LIS_constantsMod, only : LIS_CONST_PATH_LEN

  implicit none

  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: init_AGRMET    !defines the native resolution of 
                           !the input datasets and computes
                           !interpolation weights
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  public :: agrmet_struc
!EOP  

  type, public ::  agrmet_type_dec
     real                   :: ts
     integer                :: radProcessInterval
     integer                :: pcpProcessInterval
     integer                :: readgfsInterval
     real*8                 :: radProcessAlarmTime
     real*8                 :: pcpProcessAlarmTime
     real*8                 :: gfsAlarmTime
     real*8                 :: pcpclimoAlarmTime
     real*8                 :: albAlarmTime
     real*8                 :: gfracAlarmTime
     character*100          :: agrmetdir      
     character*100          :: climodir     
     character*100          :: maskfile,maskfile2,maskfile8,maskfile16,maskfile64,maskfilell
     character*100          :: terrainfile,terrainfile8,terrainfile16,terrainfile64,terrainfilll
     character*100          :: sfcntmfile     
     character*100          :: sfcalcdir
     character*100          :: mrgpcpdir
     character*100          :: clouddir
     character*100          :: ssmidir
     character*100          :: geodir
     character*100          :: gfsdir
     character*100          :: galwemdir
     character*100          :: cdmsdir
     character*100          :: cmordir
     character*100          :: analysisdir
     character*100          :: retroFileRoot
     character*20           :: first_guess_source
     integer                :: use_timestamp
     integer                :: gfs_timestamp
     real*8                 :: agrmettime1, agrmettime2
     real*8                 :: agrmetpcptime1, agrmetpcptime2
     real*8                 :: cmortime
     integer                :: pcpobswch
     integer :: pcpobsfmt ! EMK...File format for precip obs
     integer :: sfcobsfmt ! EMK...File format for sfcobs
     integer                :: pwswch
     integer                :: raswch
     integer                :: cdfs2swch
     integer                :: clswch
     integer                :: geoswch
     integer                :: cmorswch
     integer                :: imerg_swch
     integer                :: cmorminthresh
     integer                :: cmormaxthresh
     integer                :: geominthresh
     integer                :: geomaxthresh
     integer                :: imerg_t_thresh
     integer                :: gfsprecswch
     integer                :: galwemprecswch
     integer                :: razero
     real                   :: salp
     real                   :: clmult
     real                   :: mnpr
     real                   :: mxpr
     real                   :: mnpd
     real                   :: mxpd
     real                   :: cldth
     real                   :: mdmv
     real                   :: mdpe
     real                   :: ovpe
     integer                :: diff_grid,global_or_hemi 
     integer                :: cdfs2int
     real                   :: pcap 
     real                   :: wndwgt
     real                   :: minwnd
     integer                :: gridspan
     integer                :: mo1
     integer                :: mo2
     integer                :: hemi_nc(2)
     integer                :: hemi_nr(2)
     integer                :: shemi
     integer                :: nhemi
     integer                :: mi,mi2,mi3,micmor
     integer                :: imax,imax2,imax3,imaxnative,imaxsmi,imaxgp,imaxcmor
     integer                :: jmax,jmax2,jmax3,jmaxnative,jmaxsmi,jmaxgp,jmaxcmor
     real                   :: cmorminlat, cmormaxlat, cmorminlon, cmormaxlon, cmordx, cmordy
     integer                :: interp
     character*50           :: met_interp
     character*50           :: fg_gfs_interp
     character*50           :: fg_galwem_interp
     integer,allocatable    :: land(:,:,:),land2(:,:,:)
     real,allocatable       :: alt(:,:,:)
     real, allocatable      :: irad(:,:)  !sfcalc - cntnm
     real, allocatable      :: radius(:,:)
     real, allocatable      :: rlat1_nh(:),rlat1_nh2(:),rlat1_nh3(:)
     real, allocatable      :: rlat1_sh(:),rlat1_sh2(:),rlat1_sh3(:)
     real, allocatable      :: rlon1_nh(:),rlon1_nh2(:),rlon1_nh3(:)
     real, allocatable      :: rlon1_sh(:),rlon1_sh2(:),rlon1_sh3(:)
!<kluge -- cod testing>
character(len=19)      :: compute_radiation
real, allocatable      :: rlat1_nh4(:)
real, allocatable      :: rlat1_sh4(:)
real, allocatable      :: rlon1_nh4(:)
real, allocatable      :: rlon1_sh4(:)
integer, allocatable   :: n111_nh4(:)
integer, allocatable   :: n111_sh4(:)
integer, allocatable   :: n121_nh4(:)
integer, allocatable   :: n121_sh4(:)
integer, allocatable   :: n211_nh4(:)
integer, allocatable   :: n211_sh4(:)
integer, allocatable   :: n221_nh4(:)
integer, allocatable   :: n221_sh4(:)
real, allocatable      :: w111_nh4(:)
real, allocatable      :: w121_nh4(:)
real, allocatable      :: w111_sh4(:)
real, allocatable      :: w121_sh4(:)
real, allocatable      :: w211_nh4(:)
real, allocatable      :: w221_nh4(:)
real, allocatable      :: w211_sh4(:)
real, allocatable      :: w221_sh4(:)
real, allocatable      :: rlat2_nh4(:)
real, allocatable      :: rlon2_nh4(:)
real, allocatable      :: rlat2_sh4(:)
real, allocatable      :: rlon2_sh4(:)
integer, allocatable   :: n112_nh4(:)
integer, allocatable   :: n112_sh4(:)
!</kluge -- cod testing>
     integer, allocatable   :: n111_nh(:),n111_nh2(:),n111_nh3(:)
     integer, allocatable   :: n111_sh(:),n111_sh2(:),n111_sh3(:)
     integer, allocatable   :: n121_nh(:),n121_nh2(:),n121_nh3(:)
     integer, allocatable   :: n121_sh(:),n121_sh2(:),n121_sh3(:)
     integer, allocatable   :: n211_nh(:),n211_nh2(:),n211_nh3(:)
     integer, allocatable   :: n211_sh(:),n211_sh2(:),n211_sh3(:)
     integer, allocatable   :: n221_nh(:),n221_nh2(:),n221_nh3(:)
     integer, allocatable   :: n221_sh(:),n221_sh2(:),n221_sh3(:)
     integer, allocatable   :: n112cmor(:,:)
     integer, allocatable   :: n122cmor(:,:)
     integer, allocatable   :: n212cmor(:,:)
     integer, allocatable   :: n222cmor(:,:)
     real, allocatable      :: w112cmor(:,:),w122cmor(:,:)
     real, allocatable      :: w212cmor(:,:),w222cmor(:,:)
     real, allocatable      :: w111_nh(:),w111_nh2(:),w111_nh3(:)
     real, allocatable      :: w121_nh(:),w121_nh2(:),w121_nh3(:)
     real, allocatable      :: w111_sh(:),w111_sh2(:),w111_sh3(:)
     real, allocatable      :: w121_sh(:),w121_sh2(:),w121_sh3(:)
     real, allocatable      :: w211_nh(:),w211_nh2(:),w211_nh3(:)
     real, allocatable      :: w221_nh(:),w221_nh2(:),w221_nh3(:)
     real, allocatable      :: w211_sh(:),w211_sh2(:),w211_sh3(:)
     real, allocatable      :: w221_sh(:),w221_sh2(:),w221_sh3(:)
     real, allocatable      :: rlat2_nh(:),rlat2_nh2(:),rlat2_nh3(:)
     real, allocatable      :: rlon2_nh(:),rlon2_nh2(:),rlon2_nh3(:)
     real, allocatable      :: rlat2_sh(:),rlat2_sh2(:),rlat2_sh3(:)
     real, allocatable      :: rlon2_sh(:),rlon2_sh2(:),rlon2_sh3(:)
     integer, allocatable   :: n112_nh(:),n112_nh2(:),n112_nh3(:)
     integer, allocatable   :: n112_sh(:),n112_sh2(:),n112_sh3(:)
     integer,allocatable    :: smask1(:,:)
     integer, allocatable   :: smask2(:,:)
     logical                :: fillflag1
     logical                :: fillflag2
     real, allocatable      :: cliprc1(:,:)
     real, allocatable      :: clippd1(:,:)
     real, allocatable      :: clirtn1(:,:)
     real, allocatable      :: cliprc2(:,:)
     real, allocatable      :: clippd2(:,:)
     real, allocatable      :: clirtn2(:,:)
     real, allocatable      :: cliprc(:,:)
     real, allocatable      :: clippd(:,:)
     real, allocatable      :: clirtn(:,:)
     real, allocatable      :: mrgp(:,:,:)
     logical                :: pcp_start
     logical                :: pcp_ready
     character(len=6)       :: agr_bgrd_src_c
     real, allocatable      :: agr_hgt_c  ( : , : , : )
     real, allocatable      :: agr_rh_c   ( : , : , : )
     real, allocatable      :: agr_tmp_c  ( : , : , : )
     real, allocatable      :: agr_hgt_sfc_c  ( : , : )
     real, allocatable      :: agr_rh_sfc_c   ( : , : )
     real, allocatable      :: agr_tmp_sfc_c  ( : , : )
     real, allocatable      :: agr_wspd_c ( : , : )
     real, allocatable      :: agr_pres_c ( : , : )
     character(len=6)       :: agr_bgrd_src_p
     real, allocatable      :: agr_hgt_p  ( : , : , : )
     real, allocatable      :: agr_rh_p   ( : , : , : )
     real, allocatable      :: agr_tmp_p  ( : , : , : )
     real, allocatable      :: agr_hgt_sfc_p  ( : , : )
     real, allocatable      :: agr_rh_sfc_p   ( : , : )
     real, allocatable      :: agr_tmp_sfc_p  ( : , : )
     real, allocatable      :: agr_wspd_p ( : , : )
     real, allocatable      :: agr_pres_p ( : , : )

!<rm -- jim merge>
#if 0
! EMK BEGIN
     real, allocatable      :: agr_hgt_c_glb  ( : , : , : )
     real, allocatable      :: agr_rh_c_glb   ( : , : , : )
     real, allocatable      :: agr_tmp_c_glb  ( : , : , : )
     real, allocatable      :: agr_hgt_sfc_c_glb  ( : , : )
     real, allocatable      :: agr_rh_sfc_c_glb   ( : , : )
     real, allocatable      :: agr_tmp_sfc_c_glb  ( : , : )
     real, allocatable      :: agr_wspd_c_glb ( : , : )
     real, allocatable      :: agr_pres_c_glb ( : , : )

     real, allocatable      :: agr_hgt_p_glb  ( : , : , : )
     real, allocatable      :: agr_rh_p_glb   ( : , : , : )
     real, allocatable      :: agr_tmp_p_glb  ( : , : , : )
     real, allocatable      :: agr_hgt_sfc_p_glb  ( : , : )
     real, allocatable      :: agr_rh_sfc_p_glb   ( : , : )
     real, allocatable      :: agr_tmp_sfc_p_glb  ( : , : )
     real, allocatable      :: agr_wspd_p_glb ( : , : )
     real, allocatable      :: agr_pres_p_glb ( : , : )

! EMK END
#endif
!</rm -- jim merge>

     integer                :: step 
     real, allocatable      :: sfcprs(:,:,:)
     real, allocatable      :: sfcrlh(:,:,:)
     real, allocatable      :: sfcspd(:,:,:)
     real, allocatable      :: sfctmp(:,:,:)
     real, allocatable      :: lasprs(:,:)
     real, allocatable      :: lasrlh(:,:)
     real, allocatable      :: lasspd(:,:)
     real, allocatable      :: lastmp(:,:)
! EMK BEGIN
!<rm -- jim merge>
#if 0
     real, allocatable      :: sfcprs_glb(:,:,:)
#endif
!</rm -- jim merge>
     real, allocatable      :: sfcrlh_glb(:,:,:)
     real, allocatable      :: sfcspd_glb(:,:,:)
     real, allocatable      :: sfctmp_glb(:,:,:)
!<rm -- jim merge>
#if 0
     real, allocatable      :: lasprs_glb(:,:)
     real, allocatable      :: lasrlh_glb(:,:)
     real, allocatable      :: lasspd_glb(:,:)
     real, allocatable      :: lastmp_glb(:,:)

#endif
!</rm -- jim merge>
! EMK END
     integer                :: lastSfcalcHour
     integer                :: lastPcpHour
!new 
     integer            :: ncol, nrow
     integer, allocatable   :: n11_1_gfs(:)
     integer, allocatable   :: n12_1_gfs(:)
     integer, allocatable   :: n21_1_gfs(:)
     integer, allocatable   :: n22_1_gfs(:)
     integer, allocatable   :: n11_2_gfs(:,:)
     integer, allocatable   :: n12_2_gfs(:,:)
     integer, allocatable   :: n21_2_gfs(:,:)
     integer, allocatable   :: n22_2_gfs(:,:)
     integer, allocatable   :: n11_1_galwem(:)
     integer, allocatable   :: n12_1_galwem(:)
     integer, allocatable   :: n21_1_galwem(:)
     integer, allocatable   :: n22_1_galwem(:)
     integer, allocatable   :: n11_2_galwem(:,:)
     integer, allocatable   :: n12_2_galwem(:,:)
     integer, allocatable   :: n21_2_galwem(:,:)
     integer, allocatable   :: n22_2_galwem(:,:)
     integer, allocatable   :: n11_anl(:)
     integer, allocatable   :: n12_anl(:)
     integer, allocatable   :: n21_anl(:)
     integer, allocatable   :: n22_anl(:)
     integer :: mi111 ! EMK
     integer, allocatable   :: n111_anl(:) ! EMK
     real, allocatable      :: w11_1_gfs(:)
     real, allocatable      :: w12_1_gfs(:)
     real, allocatable      :: w21_1_gfs(:)
     real, allocatable      :: w22_1_gfs(:)
     real, allocatable      :: w11_2_gfs(:,:)
     real, allocatable      :: w12_2_gfs(:,:)
     real, allocatable      :: w21_2_gfs(:,:)
     real, allocatable      :: w22_2_gfs(:,:)
     real, allocatable      :: w11_1_galwem(:)
     real, allocatable      :: w12_1_galwem(:)
     real, allocatable      :: w21_1_galwem(:)
     real, allocatable      :: w22_1_galwem(:)
     real, allocatable      :: w11_2_galwem(:,:)
     real, allocatable      :: w12_2_galwem(:,:)
     real, allocatable      :: w21_2_galwem(:,:)
     real, allocatable      :: w22_2_galwem(:,:)
     real, allocatable      :: w11_anl(:)
     real, allocatable      :: w12_anl(:)
     real, allocatable      :: w21_anl(:)
     real, allocatable      :: w22_anl(:)
     integer            :: max_sfcobs
     integer            :: max_pcpobs
     integer            :: findtime1, findtime2
     real, allocatable :: metdata1(:,:)
     real, allocatable :: metdata2(:,:)

     ! EMK...IMERG settings
     character*100          :: imerg_dir
     character*10           :: imerg_product
     character*10           :: imerg_version
     integer*2              :: imerg_plp_thresh

     ! EMK NEW...For input forcing data
     character(len=50) :: forcingMapProj
     real :: forcingLowerLeftLat
     real :: forcingLowerLeftLon
     real :: forcingUpperRightLat
     real :: forcingUpperRightLon
     real :: forcingResDx
     real :: forcingResDy

     ! EMK NEW...For Bratseth scheme using GALWEM background
     real :: galwem_precip_back_err_scale_length
     real :: galwem_precip_back_sigma_b_sqr
     real :: galwem_precip_gauge_sigma_o_sqr
     real :: galwem_precip_geoprecip_err_scale_length
     real :: galwem_precip_geoprecip_sigma_o_sqr
     real :: galwem_precip_ssmi_err_scale_length
     real :: galwem_precip_ssmi_sigma_o_sqr
     real :: galwem_precip_cmorph_err_scale_length
     real :: galwem_precip_cmorph_sigma_o_sqr
     real :: galwem_precip_imerg_err_scale_length
     real :: galwem_precip_imerg_sigma_o_sqr
     real :: galwem_precip_max_dist

     real :: galwem_t2m_back_err_scale_length
     real :: galwem_t2m_back_sigma_b_sqr
     real :: galwem_t2m_stn_sigma_o_sqr
     real :: galwem_t2m_max_dist

     real :: galwem_rh2m_back_err_scale_length
     real :: galwem_rh2m_back_sigma_b_sqr
     real :: galwem_rh2m_stn_sigma_o_sqr
     real :: galwem_rh2m_max_dist

     real :: galwem_spd10m_back_err_scale_length
     real :: galwem_spd10m_back_sigma_b_sqr
     real :: galwem_spd10m_stn_sigma_o_sqr
     real :: galwem_spd10m_max_dist

     ! EMK NEW...For Bratseth scheme using GFS background
     real :: gfs_precip_back_err_scale_length
     real :: gfs_precip_back_sigma_b_sqr
     real :: gfs_precip_gauge_sigma_o_sqr
     real :: gfs_precip_geoprecip_err_scale_length
     real :: gfs_precip_geoprecip_sigma_o_sqr
     real :: gfs_precip_ssmi_err_scale_length
     real :: gfs_precip_ssmi_sigma_o_sqr
     real :: gfs_precip_cmorph_err_scale_length
     real :: gfs_precip_cmorph_sigma_o_sqr
     real :: gfs_precip_imerg_err_scale_length
     real :: gfs_precip_imerg_sigma_o_sqr
     real :: gfs_precip_max_dist

     real :: gfs_t2m_back_err_scale_length
     real :: gfs_t2m_back_sigma_b_sqr
     real :: gfs_t2m_stn_sigma_o_sqr
     real :: gfs_t2m_max_dist

     real :: gfs_rh2m_back_err_scale_length
     real :: gfs_rh2m_back_sigma_b_sqr
     real :: gfs_rh2m_stn_sigma_o_sqr
     real :: gfs_rh2m_max_dist

     real :: gfs_spd10m_back_err_scale_length
     real :: gfs_spd10m_back_sigma_b_sqr
     real :: gfs_spd10m_stn_sigma_o_sqr
     real :: gfs_spd10m_max_dist

     ! EMK NEW...For Bratseth scheme.  Set dynamically based on available
     ! background field
     real :: bratseth_precip_back_err_scale_length
     real :: bratseth_precip_back_sigma_b_sqr
     real :: bratseth_precip_gauge_sigma_o_sqr
     real :: bratseth_precip_geoprecip_err_scale_length
     real :: bratseth_precip_geoprecip_sigma_o_sqr
     real :: bratseth_precip_ssmi_err_scale_length
     real :: bratseth_precip_ssmi_sigma_o_sqr
     real :: bratseth_precip_cmorph_err_scale_length
     real :: bratseth_precip_cmorph_sigma_o_sqr
     real :: bratseth_precip_imerg_err_scale_length
     real :: bratseth_precip_imerg_sigma_o_sqr
     real :: bratseth_precip_max_dist

     real :: bratseth_t2m_back_err_scale_length
     real :: bratseth_t2m_back_sigma_b_sqr
     real :: bratseth_t2m_stn_sigma_o_sqr
     real :: bratseth_t2m_max_dist

     real :: bratseth_rh2m_back_err_scale_length
     real :: bratseth_rh2m_back_sigma_b_sqr
     real :: bratseth_rh2m_stn_sigma_o_sqr
     real :: bratseth_rh2m_max_dist

     real :: bratseth_spd10m_back_err_scale_length
     real :: bratseth_spd10m_back_sigma_b_sqr
     real :: bratseth_spd10m_stn_sigma_o_sqr
     real :: bratseth_spd10m_max_dist

     integer :: galwem_res ! EMK for GALWEM 10-km     
     integer :: gfs_filename_version ! EMK for new GFS filename convention

     ! EMK Add OBA option
     integer :: oba_switch
     integer :: skip_backqc
     integer :: skip_superstatqc

     ! EMK Add WWMCA GRIB1 option
     integer :: read_wwmca_grib1

     ! EMK Add GFS-to-GALWEM bias correction
     integer :: back_bias_corr
     real, allocatable :: pcp_back_bias_ratio(:,:)
     ! EMK Add NRT bias correction toward IMERG-Final Run
     ! (back_bias_corr == 2)
     character(LIS_CONST_PATH_LEN) :: gfs_nrt_bias_ratio_file
     character(LIS_CONST_PATH_LEN) :: galwem_nrt_bias_ratio_file
     real, allocatable :: gfs_nrt_bias_ratio(:,:)
     real, allocatable :: galwem_nrt_bias_ratio(:,:)
     integer :: pcp_back_bias_ratio_month
  end type agrmet_type_dec

  type(agrmet_type_dec), allocatable :: agrmet_struc(:)
contains
!BOP
!
! !ROUTINE: init_AGRMET
! \label{init_AGRMET}
!
! !REVISION HISTORY: 
! 25Jul2005: Sujay Kumar; Initial Specification
! 14Aug2008: Chris Franks; Add pcpblacklistfile item
! 31 MAR 2010 Add handling of various precip forcing
!             resolutions FOR LATLON LIS GRID ONLY....
!             ................see in code comments
!             for where..........Michael Shaw/WXE
! 30 AUG 2010 Turned "on" handling of non-ECE
!             (i.e. non- equatorial cylindrical, so
!             lambert, polar stereo) LIS domains
! 10 SEP 2010 Added variables for maximum number of
!             sfc & precip obs..Chris Franks/16WS/WXE/SEMS
! 20Jun2011: Ryan Ruhge; Add CMORPH
! 
! !INTERFACE:
  subroutine init_AGRMET(findex)

! !USES: 
    use LIS_coreMod,    only : LIS_rc
    use LIS_timeMgrMod, only : LIS_update_timestep
    use LIS_logMod,     only : LIS_logunit, LIS_endrun
    use LIS_snowMod,    only : LIS_snow_setup
    use LIS_logMod,     only : LIS_abort, LIS_endrun, &
           LIS_getNextUnitNumber, LIS_releaseUnitNumber
    use LIS_spatialDownscalingMod, only : LIS_init_pcpclimo_native
    use LIS_pluginIndices

    implicit none
! !ARGUMENTS: 
    integer, intent(in) :: findex
! 
! !DESCRIPTION: 
!  Defines the native resolution of the input forcing for AFWA
!  data. The grid description arrays are based on the decoding
!  schemes used by NCEP and followed in the LIS interpolation
!  schemes \ref{interp}. The routine first reads the runtime 
!  configurable options specific to various AGRMET analyses. The 
!  interpolation weights to perform the geospatial transforms 
!  to the LIS running domain is computed to be used later. This
!  routine also invokes calls to read all the static data on the
!  AGRMET grid such as landmask, terrain, and the spreading 
!  radii used for the barnes analyses. All the time-based alarms
!  to read various AGRMET data is also initialize in this routine. 
!
!  The routines invoked are: 
!  \begin{description}
!   \item[readcrd\_agrmet](\ref{readcrd_agrmet}) \newline
!     reads the runtime options specified for AGRMET algorithms
!   \item[AGRMET\_readmask](\ref{AGRMET_readmask}) \newline
!     reads the AGRMET landmask
!   \item[AGRMET\_readterrain](\ref{AGRMET_readterrain}) \newline
!     reads the AGRMET terrain
!   \item[AGRMET\_readpcpcntm](\ref{AGRMET_readpcpcntm}) \newline
!     reads the spreading radii for the precipitation analysis
!     on the AGRMET grid
!   \item[AGRMET\_read\_sfcalccntm](\ref{AGRMET_read_sfcalccntm}) \newline
!     reads the spreading radii for the surface calculations
!     on the AGRMET grid
!   \item[bilinear\_interp\_input](\ref{bilinear_interp_input}) \newline
!    computes the neighbor and  weights for bilinear interpolation
!   \item[neighbor\_interp\_input](\ref{neighbor_interp_input}) \newline
!    computes the neighbor, weights for neighbor interpolation
!   \item[conserv\_interp\_input](\ref{conserv_interp_input}) \newline
!    computes the neighbor, weights for conservative interpolation
!   \item[polarToLatLon](\ref{polarToLatLon}) \newline
!    computes lat lon values for AGRMET grid points
!   \item[AGRMET\_read\_pcpclimodata](\ref{AGRMET_read_pcpclimodata}) \newline
!    reads the precip climatology files
!  \end{description}
!EOP

    real :: gridDesci(50),gridDesci2(50),gridDesci3(50)
!<kluge -- cod testing>
real :: griddesci4(50)
!</kluge -- cod testing>
    real :: gridDesco(50)
    real :: gridDesci_glb(50)
    real :: gridDesciCmor(50)    
    real :: xmeshl,xmeshl2,xmeshl3
    real :: xpnmcaf,xpnmcaf2,xpnmcaf3
    real :: ypnmcaf,ypnmcaf2,ypnmcaf3
    real :: xmesh, orient,xi1,xj1,xmesh2, orient2,xi12,xj12,xmesh3, orient3,xi13,xj13
    real :: alat1,alon1,alat12,alon12,alat13,alon13
!<kluge -- cod testing>
real :: xi14,xj14,xmesh4,orient4,alat14,alon14
!</kluge -- cod testing>
    integer :: ihemi
    integer :: n 
    integer :: kprs
    byte, allocatable :: buffer(:,:,:) 
    character*9                   :: cstat
    character*100                 :: file_name,file_nam
    character*255                 :: message(20)
    integer                       :: rec_length
    integer                       :: istat
    integer                       :: istat1
    integer                       :: ftn

    external :: readcrd_agrmet
    external :: AGRMET_readmask
    external :: AGRMET_readterrain
    external :: AGRMET_readpcpcntm
    external :: AGRMET_read_sfcalccntm
    external :: polarToLatLon
    external :: bilinear_interp_input_withgrid
    external :: neighbor_interp_input_withgrid
    external :: AGRMET_read_pcpclimodata
    external :: gfs_reset_interp_input
    external :: galwem_reset_interp_input
    external :: conserv_interp_input
    external :: bilinear_interp_input
    external :: upscaleByAveraging_input

    allocate(agrmet_struc(LIS_rc%nnest))
    call readcrd_agrmet()

    do n=1, LIS_rc%nnest
       agrmet_struc(n)%ts = 3*60*60.0 
       call LIS_update_timestep(LIS_rc, n, agrmet_struc(n)%ts)
    enddo

    LIS_rc%met_nf(findex) = 9 

    kprs = 13
    call LIS_snow_setup

    do n=1,LIS_rc%nnest
       allocate(agrmet_struc(n)%metdata1(LIS_rc%met_nf(findex),&
                LIS_rc%ngrid(n)))
       allocate(agrmet_struc(n)%metdata2(LIS_rc%met_nf(findex),&
                LIS_rc%ngrid(n)))

       agrmet_struc(n)%findtime1 = 0 
       agrmet_struc(n)%findtime2 = 0 
       if(LIS_rc%runmode.eq.LIS_agrmetrunId) then ! offline mode

! Michael Shaw - If native grid not same as geoprecip grid, will need to do some
! seperate handling of geoprecip, including reading and use of it's
! land mask (used in AGRMET_geoest.F90)
          agrmet_struc(n)%met_interp = LIS_rc%met_interp(findex)
          
          ! EMK Initialize dimensions
          agrmet_struc(n)%imax2 = 0
          agrmet_struc(n)%jmax2 = 0

          if(agrmet_struc(n)%imaxgp /= agrmet_struc(n)%imaxnative)then
             agrmet_struc(n)%imax2=agrmet_struc(n)%imaxgp
             agrmet_struc(n)%jmax2=agrmet_struc(n)%jmaxgp
             agrmet_struc(n)%mi2 = agrmet_struc(n)%imax2*agrmet_struc(n)%jmax2
             if(agrmet_struc(n)%imax2 == 512)then
                xmeshl2 = 47.625
                xpnmcaf2 = 257
                ypnmcaf2 = 257
                agrmet_struc(n)%maskfile2=agrmet_struc(n)%maskfile8
             elseif(agrmet_struc(n)%imax2 == 1024)then
                xmeshl2 = 47.625/2
                xpnmcaf2 = 513
                ypnmcaf2 = 513
                agrmet_struc(n)%maskfile2=agrmet_struc(n)%maskfile16
             elseif(agrmet_struc(n)%imax2 == 4096)then
                xmeshl2 = 47.625/8
                xpnmcaf2 = 2049  
                ypnmcaf2 = 2049
                agrmet_struc(n)%maskfile2=agrmet_struc(n)%maskfile64
             endif
             if(agrmet_struc(n)%imax2 == 3600)then
                file_name=agrmet_struc(n)%maskfilell

! Michael Shaw - Need to read this data into byte sized buffer, then copy into land structure

                allocate(buffer(agrmet_struc(n)%imax2,agrmet_struc(n)%jmax2,1), &
                  agrmet_struc(n)%land2(agrmet_struc(n)%imax2,agrmet_struc(n)%jmax2,1))
                rec_length = agrmet_struc(n)%imax2 * agrmet_struc(n)%jmax2
                ftn = LIS_getNextUnitNumber()
                open( ftn, file=trim(file_name), form='unformatted', &
                access='direct', recl=rec_length, iostat=istat )
                if( istat .eq. 0 ) then
                   read( ftn, rec=1, iostat=istat ) buffer !agrmet_struc(n)%land(:,:,1)
                   if( istat .ne. 0 ) then
                      write(cstat,'(i9)',iostat=istat1) istat
                      if( istat1 .eq. 0 )then
                         message(4) = '  status = '//trim(cstat)
                      endif
                      call LIS_abort( message )
                      call LIS_endrun
                   endif
                else
                   write(cstat,'(i9)',iostat=istat1) istat
                   if( istat1 .eq. 0 )then
                      message(1) = '  status = '//trim(cstat)
                   endif
                   call LIS_abort( message )
                   stop
                endif
                agrmet_struc(n)%land2(:,:,1)=buffer(:,:,1)
                call LIS_releaseUnitNumber(ftn)
                deallocate(buffer)
             endif
          endif

! Michael Shaw - Need these structures for any different ssmi/s dataset (i.e., one 
! that's on a grid different from the native agrmet grid fro the run 

          ! EMK Initalize variables
          agrmet_struc(n)%imax3 = 0
          agrmet_struc(n)%jmax3 = 0

          if(agrmet_struc(n)%imaxsmi /= agrmet_struc(n)%imaxnative)then
             agrmet_struc(n)%imax3=agrmet_struc(n)%imaxsmi
             agrmet_struc(n)%jmax3=agrmet_struc(n)%jmaxsmi
             agrmet_struc(n)%mi3=agrmet_struc(n)%imax3*agrmet_struc(n)%jmax3
             if(agrmet_struc(n)%imax3 == 512)then
                xmeshl3 = 47.625
                xpnmcaf3 = 257
                ypnmcaf3 = 257
             elseif(agrmet_struc(n)%imax3 == 1024)then
                xmeshl3 = 47.625/2
                xpnmcaf3 = 513
                ypnmcaf3 = 513
             elseif(agrmet_struc(n)%imax3 == 4096)then
                xmeshl3 = 47.625/8
                xpnmcaf3 = 2049 
                ypnmcaf3 = 2049 
             endif
          endif

! Michael Shaw - For future interp_agrmetvars calls, make sure the indicators of whether global input and 
! diff_grid, period, are set to 0 while the cloud and radiation stuff is interpolated - 
! this will be set accordingly in makest if working with multiple input grids
! Also, set up the "standard" basis grid stuff here

          agrmet_struc(n)%global_or_hemi=0
          agrmet_struc(n)%diff_grid=0
          agrmet_struc(n)%imax=agrmet_struc(n)%imaxnative
          agrmet_struc(n)%jmax=agrmet_struc(n)%jmaxnative
          if(agrmet_struc(n)%imax.eq.512)then ! this is for the basis grid
             agrmet_struc(n)%maskfile=agrmet_struc(n)%maskfile8
             agrmet_struc(n)%terrainfile=agrmet_struc(n)%terrainfile8
             xmeshl = 47.625
             xpnmcaf = 257
             ypnmcaf = 257
          elseif(agrmet_struc(n)%imax.eq.1024)then
             agrmet_struc(n)%maskfile=agrmet_struc(n)%maskfile16
             agrmet_struc(n)%terrainfile=agrmet_struc(n)%terrainfile16
             xmeshl = 47.625/2
             xpnmcaf = 513
             ypnmcaf = 513
          elseif(agrmet_struc(n)%imax.eq.4096)then
             agrmet_struc(n)%maskfile=agrmet_struc(n)%maskfile64
             agrmet_struc(n)%terrainfile=agrmet_struc(n)%terrainfile64
             xmeshl = 47.625/8
             xpnmcaf = 2049 
             ypnmcaf = 2049 
          endif
          agrmet_struc(n)%interp = 1
          agrmet_struc(n)%mi = agrmet_struc(n)%imax*agrmet_struc(n)%jmax
          allocate(agrmet_struc(n)%land(agrmet_struc(n)%imax,&
               agrmet_struc(n)%jmax,2))
          allocate(agrmet_struc(n)%alt(agrmet_struc(n)%imax,&
               agrmet_struc(n)%jmax,2))

! Michael Shaw - Set up all the structures that are processor dependent - agrmet stuff over LIS grid

          allocate(agrmet_struc(n)%radius(LIS_rc%lnc(n), LIS_rc%lnr(n)))
          allocate(agrmet_struc(n)%irad(LIS_rc%lnc(n), LIS_rc%lnr(n)))
          allocate(agrmet_struc(n)%agr_hgt_c(LIS_rc%lnc(n), LIS_rc%lnr(n),kprs))
          allocate(agrmet_struc(n)%agr_rh_c(LIS_rc%lnc(n), LIS_rc%lnr(n),kprs))
          allocate(agrmet_struc(n)%agr_tmp_c(LIS_rc%lnc(n), LIS_rc%lnr(n),kprs))
          allocate(agrmet_struc(n)%agr_hgt_sfc_c(LIS_rc%lnc(n), LIS_rc%lnr(n)))
          allocate(agrmet_struc(n)%agr_rh_sfc_c(LIS_rc%lnc(n), LIS_rc%lnr(n)))
          allocate(agrmet_struc(n)%agr_tmp_sfc_c(LIS_rc%lnc(n), LIS_rc%lnr(n)))
          allocate(agrmet_struc(n)%agr_wspd_c(LIS_rc%lnc(n), LIS_rc%lnr(n)))
          allocate(agrmet_struc(n)%agr_pres_c(LIS_rc%lnc(n), LIS_rc%lnr(n)))
       
          allocate(agrmet_struc(n)%agr_hgt_p(LIS_rc%lnc(n), LIS_rc%lnr(n),kprs))
          allocate(agrmet_struc(n)%agr_rh_p(LIS_rc%lnc(n), LIS_rc%lnr(n),kprs))
          allocate(agrmet_struc(n)%agr_tmp_p(LIS_rc%lnc(n), LIS_rc%lnr(n),kprs))
          allocate(agrmet_struc(n)%agr_hgt_sfc_p(LIS_rc%lnc(n), LIS_rc%lnr(n)))
          allocate(agrmet_struc(n)%agr_rh_sfc_p(LIS_rc%lnc(n), LIS_rc%lnr(n)))
          allocate(agrmet_struc(n)%agr_tmp_sfc_p(LIS_rc%lnc(n), LIS_rc%lnr(n)))
          allocate(agrmet_struc(n)%agr_wspd_p(LIS_rc%lnc(n), LIS_rc%lnr(n)))
          allocate(agrmet_struc(n)%agr_pres_p(LIS_rc%lnc(n), LIS_rc%lnr(n)))
          agrmet_struc(n)%agr_bgrd_src_c = "NULL"
          agrmet_struc(n)%agr_hgt_c = 0.0
          agrmet_struc(n)%agr_rh_c = 0.0
          agrmet_struc(n)%agr_tmp_c = 0.0
          agrmet_struc(n)%agr_wspd_c = 0.0       
          agrmet_struc(n)%agr_pres_c = 0.0       

          agrmet_struc(n)%agr_bgrd_src_p = "NULL"
          agrmet_struc(n)%agr_hgt_p = 0.0
          agrmet_struc(n)%agr_rh_p = 0.0
          agrmet_struc(n)%agr_tmp_p = 0.0
          agrmet_struc(n)%agr_wspd_p = 0.0
          agrmet_struc(n)%agr_pres_p = 0.0

!<rm -- jim merge>
#if 0
! EMK BEGIN
          allocate(agrmet_struc(n)%agr_hgt_c_glb(LIS_rc%gnc(n), LIS_rc%gnr(n),kprs))
          allocate(agrmet_struc(n)%agr_rh_c_glb(LIS_rc%gnc(n), LIS_rc%gnr(n),kprs))
          allocate(agrmet_struc(n)%agr_tmp_c_glb(LIS_rc%gnc(n), LIS_rc%gnr(n),kprs))
          allocate(agrmet_struc(n)%agr_hgt_sfc_c_glb(LIS_rc%gnc(n), LIS_rc%gnr(n)))
          allocate(agrmet_struc(n)%agr_rh_sfc_c_glb(LIS_rc%gnc(n), LIS_rc%gnr(n)))
          allocate(agrmet_struc(n)%agr_tmp_sfc_c_glb(LIS_rc%gnc(n), LIS_rc%gnr(n)))

          allocate(agrmet_struc(n)%agr_wspd_c_glb(LIS_rc%gnc(n), LIS_rc%gnr(n)))
          allocate(agrmet_struc(n)%agr_pres_c_glb(LIS_rc%gnc(n), LIS_rc%gnr(n)))
       
          allocate(agrmet_struc(n)%agr_hgt_p_glb(LIS_rc%gnc(n), LIS_rc%gnr(n),kprs))
          allocate(agrmet_struc(n)%agr_rh_p_glb(LIS_rc%gnc(n), LIS_rc%gnr(n),kprs))
          allocate(agrmet_struc(n)%agr_tmp_p_glb(LIS_rc%gnc(n), LIS_rc%gnr(n),kprs))
          allocate(agrmet_struc(n)%agr_hgt_sfc_p_glb(LIS_rc%gnc(n), LIS_rc%gnr(n)))
          allocate(agrmet_struc(n)%agr_rh_sfc_p_glb(LIS_rc%gnc(n), LIS_rc%gnr(n)))
          allocate(agrmet_struc(n)%agr_tmp_sfc_p_glb(LIS_rc%gnc(n), LIS_rc%gnr(n)))
          allocate(agrmet_struc(n)%agr_wspd_p_glb(LIS_rc%gnc(n), LIS_rc%gnr(n)))
          allocate(agrmet_struc(n)%agr_pres_p_glb(LIS_rc%gnc(n), LIS_rc%gnr(n)))
          agrmet_struc(n)%agr_hgt_c_glb = 0.0
          agrmet_struc(n)%agr_rh_c_glb = 0.0
          agrmet_struc(n)%agr_tmp_c_glb = 0.0
          agrmet_struc(n)%agr_hgt_sfc_c_glb = 0.0
          agrmet_struc(n)%agr_rh_sfc_c_glb = 0.0
          agrmet_struc(n)%agr_tmp_sfc_c_glb = 0.0
          agrmet_struc(n)%agr_wspd_c_glb = 0.0       
          agrmet_struc(n)%agr_pres_c_glb = 0.0       

          agrmet_struc(n)%agr_hgt_p_glb = 0.0
          agrmet_struc(n)%agr_rh_p_glb = 0.0
          agrmet_struc(n)%agr_tmp_p_glb = 0.0
          agrmet_struc(n)%agr_hgt_sfc_p_glb = 0.0
          agrmet_struc(n)%agr_rh_sfc_p_glb = 0.0
          agrmet_struc(n)%agr_tmp_sfc_p_glb = 0.0
          agrmet_struc(n)%agr_wspd_p_glb = 0.0
          agrmet_struc(n)%agr_pres_p_glb = 0.0
! EMK END
#endif
!</rm -- jim merge>
        
! Michael Shaw - Use standard routines for standard grid stuff - will do the non-standard landmask
! stuff in the present definition routine 
         
          call AGRMET_readmask(n)
          call AGRMET_readterrain(n)
          call AGRMET_readpcpcntm(n)
          call AGRMET_read_sfcalccntm(n)
       
! Michael Shaw - The grid span (whether it spans across the equator or is limited
! to either hemisphere) is determined. This information is currently 
! determined if for a run domain defined in the lat/lon mode. 
! Currently we can run on a regular LIS lat/lon, the AFWA lat/lon, polar stereo,
! or lambert conformal.

          if(LIS_rc%gridDesc(n,1) .eq.0) then !latlon domain
             allocate(agrmet_struc(n)%smask1(LIS_rc%lnc(n),LIS_rc%lnr(n)))
             allocate(agrmet_struc(n)%smask2(LIS_rc%lnc(n),LIS_rc%lnr(n)))
             agrmet_struc(n)%fillflag1 = .true. 
             agrmet_struc(n)%fillflag2 = .true.
             if(LIS_rc%gridDesc(n,4).ge.0.and.LIS_rc%gridDesc(n,7).ge.0) then 
                agrmet_struc(n)%gridspan = 1 ! NH only
                agrmet_struc(n)%shemi = 1
                agrmet_struc(n)%nhemi = 1
             elseif(LIS_rc%gridDesc(n,4).le.0.and.LIS_rc%gridDesc(n,7).le.0) then 
                agrmet_struc(n)%gridspan = 2 ! SH only 
                agrmet_struc(n)%shemi = 2
                agrmet_struc(n)%nhemi = 2
             else
                agrmet_struc(n)%gridspan = 3 ! NH and SH
                agrmet_struc(n)%shemi = 1 
                agrmet_struc(n)%nhemi = 2
             endif
          elseif(LIS_rc%gridDesc(n,1) .ne. 0) then !lambert or polar stereo domain
             allocate(agrmet_struc(n)%smask1(LIS_rc%lnc(n),LIS_rc%lnr(n)))
             allocate(agrmet_struc(n)%smask2(LIS_rc%lnc(n),LIS_rc%lnr(n)))
             agrmet_struc(n)%fillflag1 = .true.
             agrmet_struc(n)%fillflag2 = .true.
             if(LIS_rc%gridDesc(n,4).ge.0) then
                agrmet_struc(n)%gridspan = 1 ! NH only
                agrmet_struc(n)%shemi = 1
                agrmet_struc(n)%nhemi = 1
             elseif(LIS_rc%gridDesc(n,4).le.0) then
                agrmet_struc(n)%gridspan = 2 ! SH only
                agrmet_struc(n)%shemi = 2
                agrmet_struc(n)%nhemi = 2
             endif
          endif

          if(LIS_rc%gridDesc(n,1).eq.0)then
             agrmet_struc(n)%hemi_nc = nint((LIS_rc%gridDesc(n,8)-LIS_rc%gridDesc(n,5))&
                  /LIS_rc%gridDesc(n,9))+1
             if(agrmet_struc(n)%gridspan.eq.1) then 
                agrmet_struc(n)%hemi_nr(1) = nint((LIS_rc%gridDesc(n,7)-&
                     LIS_rc%gridDesc(n,4))/LIS_rc%gridDesc(n,10))+1
                agrmet_struc(n)%hemi_nr(2) = 0 
                agrmet_struc(n)%mo1 = agrmet_struc(n)%hemi_nc(1)*&
                     agrmet_struc(n)%hemi_nr(1)
                agrmet_struc(n)%mo2 = 0 
             elseif(agrmet_struc(n)%gridspan.eq.2) then 
                agrmet_struc(n)%hemi_nr(1) = 0 
                agrmet_struc(n)%hemi_nr(2) = nint((LIS_rc%gridDesc(n,7)-&
                     LIS_rc%gridDesc(n,4))/LIS_rc%gridDesc(n,10)+1)
                agrmet_struc(n)%mo1 = 0 
                agrmet_struc(n)%mo2 = agrmet_struc(n)%hemi_nc(2)*&
                     agrmet_struc(n)%hemi_nr(2)
             else
                agrmet_struc(n)%hemi_nr(1) = nint((LIS_rc%gridDesc(n,7)-&
                     LIS_rc%gridDesc(n,10)/2)/LIS_rc%gridDesc(n,10)+1)
                agrmet_struc(n)%hemi_nr(2) = nint((-LIS_rc%gridDesc(n,10)/2-&
                     LIS_rc%gridDesc(n,4))/LIS_rc%gridDesc(n,10)+1) 
                agrmet_struc(n)%mo1 = agrmet_struc(n)%hemi_nc(1)*&
                     agrmet_struc(n)%hemi_nr(1)
                agrmet_struc(n)%mo2 = agrmet_struc(n)%hemi_nc(2)*&
                     agrmet_struc(n)%hemi_nr(2)
             endif
             elseif(LIS_rc%gridDesc(n,1).ne.0)then
                agrmet_struc(n)%hemi_nc = LIS_rc%gridDesc(n,2)
                if(agrmet_struc(n)%gridspan.eq.1) then
                   agrmet_struc(n)%hemi_nr(1) = LIS_rc%gridDesc(n,3)
                   agrmet_struc(n)%hemi_nr(2) = 0
                   agrmet_struc(n)%mo1 = agrmet_struc(n)%hemi_nc(1)*&
                        agrmet_struc(n)%hemi_nr(1)
                   agrmet_struc(n)%mo2 = 0
                elseif(agrmet_struc(n)%gridspan.eq.2) then
                   agrmet_struc(n)%hemi_nr(1) = 0
                   agrmet_struc(n)%hemi_nr(2) = LIS_rc%gridDesc(n,3)
                   agrmet_struc(n)%mo1 = 0
                   agrmet_struc(n)%mo2 = agrmet_struc(n)%hemi_nc(2)*&
                        agrmet_struc(n)%hemi_nr(2)
                else
                   agrmet_struc(n)%hemi_nr(1) = LIS_rc%gridDesc(n,3)
                   agrmet_struc(n)%hemi_nr(2) = LIS_rc%gridDesc(n,3)
                   agrmet_struc(n)%mo1 = agrmet_struc(n)%hemi_nc(1)*&
                        agrmet_struc(n)%hemi_nr(1)
                   agrmet_struc(n)%mo2 = agrmet_struc(n)%hemi_nc(2)*&
                        agrmet_struc(n)%hemi_nr(2)
                endif
             endif ! whether latlon or not
             
             gridDesco = 0 
             gridDesci = 0 
             do ihemi = agrmet_struc(n)%shemi,agrmet_struc(n)%nhemi

                if(LIS_rc%gridDesc(n,1).eq.0)gridDesco(1) = 0
                if(LIS_rc%gridDesc(n,1).ne.0)gridDesco(1) = LIS_rc%gridDesc(n,1)
                gridDesco(2) = agrmet_struc(n)%hemi_nc(ihemi)
                gridDesco(3) = agrmet_struc(n)%hemi_nr(ihemi)
                gridDesco(5) = LIS_rc%gridDesc(n,5)
                gridDesco(8) = LIS_rc%gridDesc(n,8)
                gridDesco(6) = LIS_rc%gridDesc(n,6)
                gridDesco(9) = LIS_rc%gridDesc(n,9)
                gridDesco(10) = LIS_rc%gridDesc(n,10)
                if(LIS_rc%gridDesc(n,1).ne.0)gridDesco(11) = LIS_rc%gridDesc(n,11)
                gridDesco(20) = 255
                if(agrmet_struc(n)%gridspan.eq.1.or.&
                     agrmet_struc(n)%gridspan.eq.2) then 
                   gridDesco(4) = LIS_rc%gridDesc(n,4)
                   gridDesco(7) = LIS_rc%gridDesc(n,7)
                else
                   if(LIS_rc%gridDesc(n,1).eq.0)then
                      if(ihemi.eq.1) then 
                         gridDesco(4) = LIS_rc%gridDesc(n,9)/2
                         gridDesco(7) = LIS_rc%gridDesc(n,7)
                      else
                         gridDesco(4) = LIS_rc%gridDesc(n,4)
                         gridDesco(7) = -LIS_rc%gridDesc(n,9)/2
                      endif
                   endif
                endif

! Michael Shaw - The input grids might vary, so handle that here if so (don't use diff_grid as condition -
! set that to 0 before...could do that a variety of ways, but stick with this for now)

! Geoprecip

                   if(agrmet_struc(n)%imaxgp /= agrmet_struc(n)%imaxnative .and. &
                     agrmet_struc(n)%imaxgp /= 3600 )then

! Michael Shaw - Added this conditional for now while deciding what to do with land masks for 64th mesh,
! or whether it makes a difference (wherther it's any different from 16th, really)

                      if(agrmet_struc(n)%imaxnative /= 1024)then 
                        if(ihemi==agrmet_struc(n)%shemi)then
                           allocate(agrmet_struc(n)%land2(agrmet_struc(n)%imax2,agrmet_struc(n)%jmax2,2))
                           rec_length=agrmet_struc(n)%imax2*agrmet_struc(n)%jmax2*4
                           file_name=agrmet_struc(n)%maskfile2
                        endif
                        ftn = LIS_getNextUnitNumber()

! Michael Shaw - ASSUMING THAT shemi=1 IS ALWAYS NORTHERN HEMISPHERE, shemi=2 SOUTHERN HEMISPHERE, HERE

                        if(ihemi==1)file_nam=trim(file_name)//'_nh'
                        if(ihemi==2)file_nam=trim(file_name)//'_sh'
                        open( ftn, file=trim(file_nam), form='unformatted', &
                          access='direct', recl=rec_length, iostat=istat )
                        if( istat .eq. 0 ) then
                           read( ftn, rec=1, iostat=istat ) agrmet_struc(n)%land2(:,:,ihemi)
                           if( istat .ne. 0 ) then
                              write(cstat,'(i9)',iostat=istat1) istat
                              if( istat1 .eq. 0 )then
                                 message(4) = '  status = '//trim(cstat)
                              endif
                              call LIS_abort( message )
                              call LIS_endrun
                           endif
                        else
                           write(cstat,'(i9)',iostat=istat1) istat
                           if( istat1 .eq. 0 )then
                              message(1) = '  status = '//trim(cstat)
                           endif
                           call LIS_abort( message )
                           stop
                        endif
                        call LIS_releaseUnitNumber(ftn)
                      endif ! Michael Shaw - the "temporary" "is this a 64th mesh" conditional
               
                      if(ihemi.eq.1) then
                         xmesh2 = xmeshl2
                         orient2 = 100.0
                      else
                         xmesh2 = -xmeshl2
                         orient2 = 280.0
                      endif
                      xj12 = float(1)-ypnmcaf2
                      xi12 = float(1)-xpnmcaf2

                      call polarToLatLon(xi12,xj12,xmesh2,orient2,alat12,alon12)
  
                      gridDesci2 = 0
                      gridDesci2(1) = 5
                      gridDesci2(2) = agrmet_struc(n)%imax2
                      gridDesci2(3) = agrmet_struc(n)%jmax2
                      gridDesci2(4) = alat12
                      gridDesci2(5) = alon12
                      gridDesci2(6) = 8
                      gridDesci2(7) = orient2
                      gridDesci2(8) = xmesh2
                      gridDesci2(9) = xmesh2
                      gridDesci2(10) = 0.0
                      !if(ihemi .eq.2) then
                      !   gridDesci2(11) = 128
                      !endif
                      gridDesci2(20)= 128
                      gridDesci2(11)= 128
                      gridDesci2(13)= 1
                   elseif(agrmet_struc(n)%imaxgp == 3600)then ! Michael Shaw - otherwise, if this is a 10th deg grid
                      gridDesci2 = 0
                      gridDesci2(1) = 0
                      gridDesci2(2) = agrmet_struc(n)%imax2
                      gridDesci2(3) = agrmet_struc(n)%jmax2
                      gridDesci2(4) = -89.95
                      gridDesci2(5) = -179.95 !0.05
                      gridDesci2(6) = 128
                      gridDesci2(7) = 89.95
                      gridDesci2(8) = 179.95 !359.95
                      gridDesci2(9) = 0.1
                      gridDesci2(10) =0.1
                      gridDesci2(20) = 128
                      gridDesci2(13)= 1
                   endif

! Michael Shaw - SSMI/S

                   if(agrmet_struc(n)%imaxsmi /= agrmet_struc(n)%imaxnative .and. &
                     agrmet_struc(n)%imaxsmi /= 1440)then
                      if(ihemi.eq.1) then
                         xmesh3 = xmeshl3
                         orient3 = 100.0
                      else
                         xmesh3 = -xmeshl3
                         orient3 = 280.0
                      endif
                      xj13 = float(1)-ypnmcaf3
                      xi13 = float(1)-xpnmcaf3

                      call polarToLatLon(xi13,xj13,xmesh3,orient3,alat13,alon13)
 
                      gridDesci3 = 0
                      gridDesci3(1) = 5
                      gridDesci3(2) = agrmet_struc(n)%imax3
                      gridDesci3(3) = agrmet_struc(n)%jmax3
                      gridDesci3(4) = alat13
                      gridDesci3(5) = alon13
                      gridDesci3(6) = 8
                      gridDesci3(7) = orient3
                      gridDesci3(8) = xmesh3
                      gridDesci3(9) = xmesh3
                      gridDesci3(10) = 0.0
                      !if(ihemi .eq.2) then
                      !   gridDesci3(11) = 128
                      !endif
                      gridDesci3(20)= 128
                      gridDesci3(11)= 128
                      gridDesci3(13) = 1
                   elseif(agrmet_struc(n)%imaxsmi == 1440)then
                      gridDesci3 = 0
                      gridDesci3(1) = 0
                      gridDesci3(2) = agrmet_struc(n)%imax3
                      gridDesci3(3) = agrmet_struc(n)%jmax3
                      gridDesci3(4) = -89.875
                      gridDesci3(5) = -179.875 !0.05
                      gridDesci3(6) = 128
                      gridDesci3(7) = 89.875
                      gridDesci3(8) = 179.875 !359.95
                      gridDesci3(9) = 0.25
                      gridDesci3(10) =0.25
                      gridDesci3(20) = 128
                      gridDesci3(13) = 1
                   endif
!<kluge -- cod testing>
 if(ihemi.eq.1) then
    xmesh4 = 47.625/2
    orient4 = 100.0
 else
    xmesh4 = -47.625/2
    orient4 = 280.0
 endif
 xj14 = float(1)-513
 xi14 = float(1)-513

 call polarToLatLon(xi14,xj14,xmesh4,orient4,alat14,alon14)

 gridDesci4 = 0
 gridDesci4(1) = 5
 gridDesci4(2) = 1024
 gridDesci4(3) = 1024
 gridDesci4(4) = alat14
 gridDesci4(5) = alon14
 gridDesci4(6) = 8
 gridDesci4(7) = orient4
 gridDesci4(8) = xmesh4
 gridDesci4(9) = xmesh4
 gridDesci4(10) = 0.0
 !if(ihemi .eq.2) then
 !   gridDesci4(11) = 128
 !endif
 gridDesci4(20)= 128
 gridDesci4(11)= 128
 gridDesci4(13) = 1
!</kluge -- cod testing>

! Do some of the orientation for the standard grid here

                if(ihemi.eq.1) then 
                   xmesh = xmeshl
                   orient = 100.0
                else
                   xmesh = -xmeshl
                   orient = 280.0
                endif
                xj1 = float(1)-ypnmcaf
                xi1 = float(1)-xpnmcaf
                
                call polarToLatLon(xi1,xj1,xmesh,orient,alat1,alon1)
                
                gridDesci = 0 
                gridDesci(1) = 5
                gridDesci(2) = agrmet_struc(n)%imax
                gridDesci(3) = agrmet_struc(n)%jmax
                gridDesci(4) = alat1
                gridDesci(5) = alon1
                gridDesci(6) = 8
                gridDesci(7) = orient
                gridDesci(8) = xmesh
                gridDesci(9) = xmesh
                gridDesci(10) = 0.0       
                !if(ihemi .eq.2) then 
                !   gridDesci(20) = 128
                !   gridDesci(11) = 128
                !else
                gridDesci(20) = 128 
                gridDesci(11) = 128 
                !endif 
                gridDesci(13) = 1  ! Michael Shaw - global grid

! Michael Shaw - Now, do the structures and get the interpolation coefficients for both the standard and nonstandard
! stuff (if it exists)                

                if(agrmet_struc(n)%gridspan.eq.1) then 
                   allocate(agrmet_struc(n)%rlat1_nh(agrmet_struc(n)%mo1))
                   allocate(agrmet_struc(n)%rlon1_nh(agrmet_struc(n)%mo1))
                   allocate(agrmet_struc(n)%n111_nh(agrmet_struc(n)%mo1))
                   allocate(agrmet_struc(n)%n121_nh(agrmet_struc(n)%mo1))
                   allocate(agrmet_struc(n)%n211_nh(agrmet_struc(n)%mo1))
                   allocate(agrmet_struc(n)%n221_nh(agrmet_struc(n)%mo1))
                   allocate(agrmet_struc(n)%w111_nh(agrmet_struc(n)%mo1))
                   allocate(agrmet_struc(n)%w121_nh(agrmet_struc(n)%mo1))
                   allocate(agrmet_struc(n)%w211_nh(agrmet_struc(n)%mo1))
                   allocate(agrmet_struc(n)%w221_nh(agrmet_struc(n)%mo1))
                   call bilinear_interp_input_withgrid(gridDesci,gridDesco,&
                        agrmet_struc(n)%mo1,agrmet_struc(n)%rlat1_nh,&
                        agrmet_struc(n)%rlon1_nh,agrmet_struc(n)%n111_nh,&
                        agrmet_struc(n)%n121_nh,agrmet_struc(n)%n211_nh,&
                        agrmet_struc(n)%n221_nh,agrmet_struc(n)%w111_nh,&
                        agrmet_struc(n)%w121_nh,agrmet_struc(n)%w211_nh,&
                        agrmet_struc(n)%w221_nh)
                   if(agrmet_struc(n)%imaxnative /= agrmet_struc(n)%imaxgp)then
                      if(agrmet_struc(n)%imax2 == 3600)then
                         gridDesci2(4)=0.05
                         gridDesci2(3)=900
                      endif
                      allocate(agrmet_struc(n)%rlat1_nh2(agrmet_struc(n)%mo1))
                      allocate(agrmet_struc(n)%rlon1_nh2(agrmet_struc(n)%mo1))
                      allocate(agrmet_struc(n)%n111_nh2(agrmet_struc(n)%mo1))
                      allocate(agrmet_struc(n)%n121_nh2(agrmet_struc(n)%mo1))
                      allocate(agrmet_struc(n)%n211_nh2(agrmet_struc(n)%mo1))
                      allocate(agrmet_struc(n)%n221_nh2(agrmet_struc(n)%mo1))
                      allocate(agrmet_struc(n)%w111_nh2(agrmet_struc(n)%mo1))
                      allocate(agrmet_struc(n)%w121_nh2(agrmet_struc(n)%mo1))
                      allocate(agrmet_struc(n)%w211_nh2(agrmet_struc(n)%mo1))
                      allocate(agrmet_struc(n)%w221_nh2(agrmet_struc(n)%mo1))
                      call bilinear_interp_input_withgrid(gridDesci2,gridDesco,&
                           agrmet_struc(n)%mo1,agrmet_struc(n)%rlat1_nh2,&
                           agrmet_struc(n)%rlon1_nh2,agrmet_struc(n)%n111_nh2,&
                           agrmet_struc(n)%n121_nh2,agrmet_struc(n)%n211_nh2,&
                           agrmet_struc(n)%n221_nh2,agrmet_struc(n)%w111_nh2,&
                           agrmet_struc(n)%w121_nh2,agrmet_struc(n)%w211_nh2,&
                           agrmet_struc(n)%w221_nh2)
                   endif ! nonstandard grid
                   allocate(agrmet_struc(n)%rlat2_nh(agrmet_struc(n)%mo1))
                   allocate(agrmet_struc(n)%rlon2_nh(agrmet_struc(n)%mo1))
                   allocate(agrmet_struc(n)%n112_nh(agrmet_struc(n)%mo1))
                   call neighbor_interp_input_withgrid(gridDesci,gridDesco,&
                        agrmet_struc(n)%mo1,agrmet_struc(n)%rlat2_nh,&
                        agrmet_struc(n)%rlon2_nh,agrmet_struc(n)%n112_nh)
                   if(agrmet_struc(n)%imaxnative /= agrmet_struc(n)%imaxgp)then
                      allocate(agrmet_struc(n)%rlat2_nh2(agrmet_struc(n)%mo1))
                      allocate(agrmet_struc(n)%rlon2_nh2(agrmet_struc(n)%mo1))
                      allocate(agrmet_struc(n)%n112_nh2(agrmet_struc(n)%mo1))
                      call neighbor_interp_input_withgrid(gridDesci2,gridDesco,&
                           agrmet_struc(n)%mo1,agrmet_struc(n)%rlat2_nh2,&
                           agrmet_struc(n)%rlon2_nh2,agrmet_struc(n)%n112_nh2)
                   endif ! nonstandard grid
                   if(agrmet_struc(n)%imaxnative /= agrmet_struc(n)%imaxsmi)then
                      if(agrmet_struc(n)%imax3 == 1440)then
                         gridDesci3(4)=0.125
                         gridDesci3(3)=360
                      endif
                      allocate(agrmet_struc(n)%rlat1_nh3(agrmet_struc(n)%mo1))
                      allocate(agrmet_struc(n)%rlon1_nh3(agrmet_struc(n)%mo1))
                      allocate(agrmet_struc(n)%n111_nh3(agrmet_struc(n)%mo1))
                      allocate(agrmet_struc(n)%n121_nh3(agrmet_struc(n)%mo1))
                      allocate(agrmet_struc(n)%n211_nh3(agrmet_struc(n)%mo1))
                      allocate(agrmet_struc(n)%n221_nh3(agrmet_struc(n)%mo1))
                      allocate(agrmet_struc(n)%w111_nh3(agrmet_struc(n)%mo1))
                      allocate(agrmet_struc(n)%w121_nh3(agrmet_struc(n)%mo1))
                      allocate(agrmet_struc(n)%w211_nh3(agrmet_struc(n)%mo1))
                      allocate(agrmet_struc(n)%w221_nh3(agrmet_struc(n)%mo1))
                      call bilinear_interp_input_withgrid(gridDesci3,gridDesco,&
                           agrmet_struc(n)%mo1,agrmet_struc(n)%rlat1_nh3,&
                           agrmet_struc(n)%rlon1_nh3,agrmet_struc(n)%n111_nh3,&
                           agrmet_struc(n)%n121_nh3,agrmet_struc(n)%n211_nh3,&
                           agrmet_struc(n)%n221_nh3,agrmet_struc(n)%w111_nh3,&
                           agrmet_struc(n)%w121_nh3,agrmet_struc(n)%w211_nh3,&
                           agrmet_struc(n)%w221_nh3)
                      allocate(agrmet_struc(n)%rlat2_nh3(agrmet_struc(n)%mo1))
                      allocate(agrmet_struc(n)%rlon2_nh3(agrmet_struc(n)%mo1))
                      allocate(agrmet_struc(n)%n112_nh3(agrmet_struc(n)%mo1))
                      call neighbor_interp_input_withgrid(gridDesci3,gridDesco,&
                           agrmet_struc(n)%mo1,agrmet_struc(n)%rlat2_nh3,&
                           agrmet_struc(n)%rlon2_nh3,agrmet_struc(n)%n112_nh3)
                   endif ! nonstandard grid
!<kluge -- cod testing>
                   allocate(agrmet_struc(n)%rlat1_nh4(agrmet_struc(n)%mo1))
                   allocate(agrmet_struc(n)%rlon1_nh4(agrmet_struc(n)%mo1))
                   allocate(agrmet_struc(n)%n111_nh4(agrmet_struc(n)%mo1))
                   allocate(agrmet_struc(n)%n121_nh4(agrmet_struc(n)%mo1))
                   allocate(agrmet_struc(n)%n211_nh4(agrmet_struc(n)%mo1))
                   allocate(agrmet_struc(n)%n221_nh4(agrmet_struc(n)%mo1))
                   allocate(agrmet_struc(n)%w111_nh4(agrmet_struc(n)%mo1))
                   allocate(agrmet_struc(n)%w121_nh4(agrmet_struc(n)%mo1))
                   allocate(agrmet_struc(n)%w211_nh4(agrmet_struc(n)%mo1))
                   allocate(agrmet_struc(n)%w221_nh4(agrmet_struc(n)%mo1))
                   call bilinear_interp_input_withgrid(gridDesci4,gridDesco,&
                        agrmet_struc(n)%mo1,agrmet_struc(n)%rlat1_nh4,&
                        agrmet_struc(n)%rlon1_nh4,agrmet_struc(n)%n111_nh4,&
                        agrmet_struc(n)%n121_nh4,agrmet_struc(n)%n211_nh4,&
                        agrmet_struc(n)%n221_nh4,agrmet_struc(n)%w111_nh4,&
                        agrmet_struc(n)%w121_nh4,agrmet_struc(n)%w211_nh4,&
                        agrmet_struc(n)%w221_nh4)
                   allocate(agrmet_struc(n)%rlat2_nh4(agrmet_struc(n)%mo1))
                   allocate(agrmet_struc(n)%rlon2_nh4(agrmet_struc(n)%mo1))
                   allocate(agrmet_struc(n)%n112_nh4(agrmet_struc(n)%mo1))
                   call neighbor_interp_input_withgrid(gridDesci4,gridDesco,&
                        agrmet_struc(n)%mo1,agrmet_struc(n)%rlat2_nh4,&
                        agrmet_struc(n)%rlon2_nh4,agrmet_struc(n)%n112_nh4)
!</kluge -- cod testing>
                elseif(agrmet_struc(n)%gridspan.eq.2) then 
                   allocate(agrmet_struc(n)%rlat1_sh(agrmet_struc(n)%mo2))
                   allocate(agrmet_struc(n)%rlon1_sh(agrmet_struc(n)%mo2))
                   allocate(agrmet_struc(n)%n111_sh(agrmet_struc(n)%mo2))
                   allocate(agrmet_struc(n)%n121_sh(agrmet_struc(n)%mo2))
                   allocate(agrmet_struc(n)%n211_sh(agrmet_struc(n)%mo2))
                   allocate(agrmet_struc(n)%n221_sh(agrmet_struc(n)%mo2))
                   allocate(agrmet_struc(n)%w111_sh(agrmet_struc(n)%mo2))
                   allocate(agrmet_struc(n)%w121_sh(agrmet_struc(n)%mo2))
                   allocate(agrmet_struc(n)%w211_sh(agrmet_struc(n)%mo2))
                   allocate(agrmet_struc(n)%w221_sh(agrmet_struc(n)%mo2))
                   call bilinear_interp_input_withgrid(gridDesci,gridDesco,&
                        agrmet_struc(n)%mo2,agrmet_struc(n)%rlat1_sh,&
                        agrmet_struc(n)%rlon1_sh,&
                        agrmet_struc(n)%n111_sh,agrmet_struc(n)%n121_sh,&
                        agrmet_struc(n)%n211_sh,agrmet_struc(n)%n221_sh,&
                        agrmet_struc(n)%w111_sh,agrmet_struc(n)%w121_sh,&
                        agrmet_struc(n)%w211_sh,agrmet_struc(n)%w221_sh)
                   if(agrmet_struc(n)%imaxnative /= agrmet_struc(n)%imaxgp)then
                      if(agrmet_struc(n)%imax2 == 3600)then
                         gridDesci2(7)=-0.05
                         gridDesci2(3)=900
                      endif
                      allocate(agrmet_struc(n)%rlat1_sh2(agrmet_struc(n)%mo2))
                      allocate(agrmet_struc(n)%rlon1_sh2(agrmet_struc(n)%mo2))
                      allocate(agrmet_struc(n)%n111_sh2(agrmet_struc(n)%mo2))
                      allocate(agrmet_struc(n)%n121_sh2(agrmet_struc(n)%mo2))
                      allocate(agrmet_struc(n)%n211_sh2(agrmet_struc(n)%mo2))
                      allocate(agrmet_struc(n)%n221_sh2(agrmet_struc(n)%mo2))
                      allocate(agrmet_struc(n)%w111_sh2(agrmet_struc(n)%mo2))
                      allocate(agrmet_struc(n)%w121_sh2(agrmet_struc(n)%mo2))
                      allocate(agrmet_struc(n)%w211_sh2(agrmet_struc(n)%mo2))
                      allocate(agrmet_struc(n)%w221_sh2(agrmet_struc(n)%mo2))
                      call bilinear_interp_input_withgrid(gridDesci2,gridDesco,&
                           agrmet_struc(n)%mo2,agrmet_struc(n)%rlat1_sh2,&
                           agrmet_struc(n)%rlon1_sh2,&
                           agrmet_struc(n)%n111_sh2,agrmet_struc(n)%n121_sh2,&
                           agrmet_struc(n)%n211_sh2,agrmet_struc(n)%n221_sh2,&
                           agrmet_struc(n)%w111_sh2,agrmet_struc(n)%w121_sh2,&
                           agrmet_struc(n)%w211_sh2,agrmet_struc(n)%w221_sh2)
                   endif ! nonstandard grid
                   allocate(agrmet_struc(n)%rlat2_sh(agrmet_struc(n)%mo2))
                   allocate(agrmet_struc(n)%rlon2_sh(agrmet_struc(n)%mo2))
                   allocate(agrmet_struc(n)%n112_sh(agrmet_struc(n)%mo2))
                   call neighbor_interp_input_withgrid(gridDesci,gridDesco,&
                        agrmet_struc(n)%mo2,agrmet_struc(n)%rlat2_sh,&
                        agrmet_struc(n)%rlon2_sh,agrmet_struc(n)%n112_sh)
                   if(agrmet_struc(n)%imaxnative /= agrmet_struc(n)%imaxgp)then
                      allocate(agrmet_struc(n)%rlat2_sh2(agrmet_struc(n)%mo2))
                      allocate(agrmet_struc(n)%rlon2_sh2(agrmet_struc(n)%mo2))
                      allocate(agrmet_struc(n)%n112_sh2(agrmet_struc(n)%mo2))
                      call neighbor_interp_input_withgrid(gridDesci2,gridDesco,&
                           agrmet_struc(n)%mo2,agrmet_struc(n)%rlat2_sh2,&
                           agrmet_struc(n)%rlon2_sh2,agrmet_struc(n)%n112_sh2)
                   endif ! nonstandard grid
                   if(agrmet_struc(n)%imaxnative /= agrmet_struc(n)%imaxsmi)then
                      if(agrmet_struc(n)%imax3 == 1440)then
                         gridDesci3(7)=-0.125
                         gridDesci3(3)=360
                      endif
                      allocate(agrmet_struc(n)%rlat1_sh3(agrmet_struc(n)%mo2))
                      allocate(agrmet_struc(n)%rlon1_sh3(agrmet_struc(n)%mo2))
                      allocate(agrmet_struc(n)%n111_sh3(agrmet_struc(n)%mo2))
                      allocate(agrmet_struc(n)%n121_sh3(agrmet_struc(n)%mo2))
                      allocate(agrmet_struc(n)%n211_sh3(agrmet_struc(n)%mo2))
                      allocate(agrmet_struc(n)%n221_sh3(agrmet_struc(n)%mo2))
                      allocate(agrmet_struc(n)%w111_sh3(agrmet_struc(n)%mo2))
                      allocate(agrmet_struc(n)%w121_sh3(agrmet_struc(n)%mo2))
                      allocate(agrmet_struc(n)%w211_sh3(agrmet_struc(n)%mo2))
                      allocate(agrmet_struc(n)%w221_sh3(agrmet_struc(n)%mo2))
                      call bilinear_interp_input_withgrid(gridDesci3,gridDesco,&
                           agrmet_struc(n)%mo2,agrmet_struc(n)%rlat1_sh3,&
                           agrmet_struc(n)%rlon1_sh3,&
                           agrmet_struc(n)%n111_sh3,agrmet_struc(n)%n121_sh3,&
                           agrmet_struc(n)%n211_sh3,agrmet_struc(n)%n221_sh3,&
                           agrmet_struc(n)%w111_sh3,agrmet_struc(n)%w121_sh3,&
                           agrmet_struc(n)%w211_sh3,agrmet_struc(n)%w221_sh3)
                      allocate(agrmet_struc(n)%rlat2_sh3(agrmet_struc(n)%mo2))
                      allocate(agrmet_struc(n)%rlon2_sh3(agrmet_struc(n)%mo2))
                      allocate(agrmet_struc(n)%n112_sh3(agrmet_struc(n)%mo2))
                      call neighbor_interp_input_withgrid(gridDesci3,gridDesco,&
                           agrmet_struc(n)%mo2,agrmet_struc(n)%rlat2_sh3,&
                           agrmet_struc(n)%rlon2_sh3,agrmet_struc(n)%n112_sh3)
                   endif ! nonstandard grid
!<kluge -- cod testing>
                   allocate(agrmet_struc(n)%rlat1_sh4(agrmet_struc(n)%mo2))
                   allocate(agrmet_struc(n)%rlon1_sh4(agrmet_struc(n)%mo2))
                   allocate(agrmet_struc(n)%n111_sh4(agrmet_struc(n)%mo2))
                   allocate(agrmet_struc(n)%n121_sh4(agrmet_struc(n)%mo2))
                   allocate(agrmet_struc(n)%n211_sh4(agrmet_struc(n)%mo2))
                   allocate(agrmet_struc(n)%n221_sh4(agrmet_struc(n)%mo2))
                   allocate(agrmet_struc(n)%w111_sh4(agrmet_struc(n)%mo2))
                   allocate(agrmet_struc(n)%w121_sh4(agrmet_struc(n)%mo2))
                   allocate(agrmet_struc(n)%w211_sh4(agrmet_struc(n)%mo2))
                   allocate(agrmet_struc(n)%w221_sh4(agrmet_struc(n)%mo2))
                   call bilinear_interp_input_withgrid(gridDesci4,gridDesco,&
                        agrmet_struc(n)%mo2,agrmet_struc(n)%rlat1_sh4,&
                        agrmet_struc(n)%rlon1_sh4,agrmet_struc(n)%n111_sh4,&
                        agrmet_struc(n)%n121_sh4,agrmet_struc(n)%n211_sh4,&
                        agrmet_struc(n)%n221_sh4,agrmet_struc(n)%w111_sh4,&
                        agrmet_struc(n)%w121_sh4,agrmet_struc(n)%w211_sh4,&
                        agrmet_struc(n)%w221_sh4)
                   allocate(agrmet_struc(n)%rlat2_sh4(agrmet_struc(n)%mo2))
                   allocate(agrmet_struc(n)%rlon2_sh4(agrmet_struc(n)%mo2))
                   allocate(agrmet_struc(n)%n112_sh4(agrmet_struc(n)%mo2))
                   call neighbor_interp_input_withgrid(gridDesci4,gridDesco,&
                        agrmet_struc(n)%mo2,agrmet_struc(n)%rlat2_sh4,&
                        agrmet_struc(n)%rlon2_sh4,agrmet_struc(n)%n112_sh4)
!</kluge -- cod testing>
                else
                   if(ihemi.eq.1) then 
                      allocate(agrmet_struc(n)%rlat1_nh(agrmet_struc(n)%mo1))
                      allocate(agrmet_struc(n)%rlon1_nh(agrmet_struc(n)%mo1))
                      allocate(agrmet_struc(n)%n111_nh(agrmet_struc(n)%mo1))
                      allocate(agrmet_struc(n)%n121_nh(agrmet_struc(n)%mo1))
                      allocate(agrmet_struc(n)%n211_nh(agrmet_struc(n)%mo1))
                      allocate(agrmet_struc(n)%n221_nh(agrmet_struc(n)%mo1))
                      allocate(agrmet_struc(n)%w111_nh(agrmet_struc(n)%mo1))
                      allocate(agrmet_struc(n)%w121_nh(agrmet_struc(n)%mo1))
                      allocate(agrmet_struc(n)%w211_nh(agrmet_struc(n)%mo1))
                      allocate(agrmet_struc(n)%w221_nh(agrmet_struc(n)%mo1))
                      if(LIS_rc%gridDesc(n,1).eq.0) then 
                         call bilinear_interp_input_withgrid(&
                              gridDesci,gridDesco,&
                              agrmet_struc(n)%mo1,agrmet_struc(n)%rlat1_nh,&
                              agrmet_struc(n)%rlon1_nh,agrmet_struc(n)%n111_nh,&
                              agrmet_struc(n)%n121_nh,agrmet_struc(n)%n211_nh,&
                              agrmet_struc(n)%n221_nh,agrmet_struc(n)%w111_nh,&
                              agrmet_struc(n)%w121_nh,agrmet_struc(n)%w211_nh,&
                              agrmet_struc(n)%w221_nh)
                      endif
                      if(agrmet_struc(n)%imaxnative /= agrmet_struc(n)%imaxgp)then
                         if(agrmet_struc(n)%imax2 == 3600 )then
                            gridDesci2(4)=0.05
                            gridDesci2(3)=900
                         endif
                         allocate(agrmet_struc(n)%rlat1_nh2(agrmet_struc(n)%mo1))
                         allocate(agrmet_struc(n)%rlon1_nh2(agrmet_struc(n)%mo1))
                         allocate(agrmet_struc(n)%n111_nh2(agrmet_struc(n)%mo1))
                         allocate(agrmet_struc(n)%n121_nh2(agrmet_struc(n)%mo1))
                         allocate(agrmet_struc(n)%n211_nh2(agrmet_struc(n)%mo1))
                         allocate(agrmet_struc(n)%n221_nh2(agrmet_struc(n)%mo1))
                         allocate(agrmet_struc(n)%w111_nh2(agrmet_struc(n)%mo1))
                         allocate(agrmet_struc(n)%w121_nh2(agrmet_struc(n)%mo1))
                         allocate(agrmet_struc(n)%w211_nh2(agrmet_struc(n)%mo1))
                         allocate(agrmet_struc(n)%w221_nh2(agrmet_struc(n)%mo1))
                         if(LIS_rc%gridDesc(n,1).eq.0) then 
                            call bilinear_interp_input_withgrid(&
                                 gridDesci2,gridDesco,&
                                 agrmet_struc(n)%mo1,agrmet_struc(n)%rlat1_nh2,&
                                 agrmet_struc(n)%rlon1_nh2,agrmet_struc(n)%n111_nh2,&
                                 agrmet_struc(n)%n121_nh2,agrmet_struc(n)%n211_nh2,&
                                 agrmet_struc(n)%n221_nh2,agrmet_struc(n)%w111_nh2,&
                                 agrmet_struc(n)%w121_nh2,agrmet_struc(n)%w211_nh2,&
                                 agrmet_struc(n)%w221_nh2)
                         endif
                      endif ! nonstandard grid
                      allocate(agrmet_struc(n)%rlat2_nh(agrmet_struc(n)%mo1))
                      allocate(agrmet_struc(n)%rlon2_nh(agrmet_struc(n)%mo1))
                      allocate(agrmet_struc(n)%n112_nh(agrmet_struc(n)%mo1))
                      if(LIS_rc%gridDesc(n,1).eq.0) then 
                         call neighbor_interp_input_withgrid(gridDesci,gridDesco,&
                              agrmet_struc(n)%mo1,agrmet_struc(n)%rlat2_nh,&
                              agrmet_struc(n)%rlon2_nh,agrmet_struc(n)%n112_nh)
                      endif

                      if(agrmet_struc(n)%imaxnative /= agrmet_struc(n)%imaxgp)then
                         allocate(agrmet_struc(n)%rlat2_nh2(agrmet_struc(n)%mo1))
                         allocate(agrmet_struc(n)%rlon2_nh2(agrmet_struc(n)%mo1))
                         allocate(agrmet_struc(n)%n112_nh2(agrmet_struc(n)%mo1))
                         call neighbor_interp_input_withgrid(gridDesci2,gridDesco,&
                              agrmet_struc(n)%mo1,agrmet_struc(n)%rlat2_nh2,&
                              agrmet_struc(n)%rlon2_nh2,agrmet_struc(n)%n112_nh2)
                      endif ! nonstandard grid
                      if(agrmet_struc(n)%imaxnative /= agrmet_struc(n)%imaxsmi)then
                         if(agrmet_struc(n)%imax3 == 1440 )then
                            gridDesci3(4)=0.125
                            gridDesci3(3)=360
                         endif
                         allocate(agrmet_struc(n)%rlat1_nh3(agrmet_struc(n)%mo1))
                         allocate(agrmet_struc(n)%rlon1_nh3(agrmet_struc(n)%mo1))
                         allocate(agrmet_struc(n)%n111_nh3(agrmet_struc(n)%mo1))
                         allocate(agrmet_struc(n)%n121_nh3(agrmet_struc(n)%mo1))
                         allocate(agrmet_struc(n)%n211_nh3(agrmet_struc(n)%mo1))
                         allocate(agrmet_struc(n)%n221_nh3(agrmet_struc(n)%mo1))
                         allocate(agrmet_struc(n)%w111_nh3(agrmet_struc(n)%mo1))
                         allocate(agrmet_struc(n)%w121_nh3(agrmet_struc(n)%mo1))
                         allocate(agrmet_struc(n)%w211_nh3(agrmet_struc(n)%mo1))
                         allocate(agrmet_struc(n)%w221_nh3(agrmet_struc(n)%mo1))
                         if(LIS_rc%gridDesc(n,1).eq.0) then 
                            call bilinear_interp_input_withgrid(gridDesci3,gridDesco,&
                                 agrmet_struc(n)%mo1,agrmet_struc(n)%rlat1_nh3,&
                                 agrmet_struc(n)%rlon1_nh3,agrmet_struc(n)%n111_nh3,&
                                 agrmet_struc(n)%n121_nh3,agrmet_struc(n)%n211_nh3,&
                                 agrmet_struc(n)%n221_nh3,agrmet_struc(n)%w111_nh3,&
                                 agrmet_struc(n)%w121_nh3,agrmet_struc(n)%w211_nh3,&
                                 agrmet_struc(n)%w221_nh3)
                         endif
                         allocate(agrmet_struc(n)%rlat2_nh3(agrmet_struc(n)%mo1))
                         allocate(agrmet_struc(n)%rlon2_nh3(agrmet_struc(n)%mo1))
                         allocate(agrmet_struc(n)%n112_nh3(agrmet_struc(n)%mo1))
                         if(LIS_rc%gridDesc(n,1).eq.0) then 
                            call neighbor_interp_input_withgrid(gridDesci3,gridDesco,&
                                 agrmet_struc(n)%mo1,agrmet_struc(n)%rlat2_nh3,&
                                 agrmet_struc(n)%rlon2_nh3,agrmet_struc(n)%n112_nh3)
                         endif
                      endif ! nonstandard grid
!<kluge -- cod testing>
                   allocate(agrmet_struc(n)%rlat1_nh4(agrmet_struc(n)%mo1))
                   allocate(agrmet_struc(n)%rlon1_nh4(agrmet_struc(n)%mo1))
                   allocate(agrmet_struc(n)%n111_nh4(agrmet_struc(n)%mo1))
                   allocate(agrmet_struc(n)%n121_nh4(agrmet_struc(n)%mo1))
                   allocate(agrmet_struc(n)%n211_nh4(agrmet_struc(n)%mo1))
                   allocate(agrmet_struc(n)%n221_nh4(agrmet_struc(n)%mo1))
                   allocate(agrmet_struc(n)%w111_nh4(agrmet_struc(n)%mo1))
                   allocate(agrmet_struc(n)%w121_nh4(agrmet_struc(n)%mo1))
                   allocate(agrmet_struc(n)%w211_nh4(agrmet_struc(n)%mo1))
                   allocate(agrmet_struc(n)%w221_nh4(agrmet_struc(n)%mo1))
                   call bilinear_interp_input_withgrid(gridDesci4,gridDesco,&
                        agrmet_struc(n)%mo1,agrmet_struc(n)%rlat1_nh4,&
                        agrmet_struc(n)%rlon1_nh4,agrmet_struc(n)%n111_nh4,&
                        agrmet_struc(n)%n121_nh4,agrmet_struc(n)%n211_nh4,&
                        agrmet_struc(n)%n221_nh4,agrmet_struc(n)%w111_nh4,&
                        agrmet_struc(n)%w121_nh4,agrmet_struc(n)%w211_nh4,&
                        agrmet_struc(n)%w221_nh4)
                   allocate(agrmet_struc(n)%rlat2_nh4(agrmet_struc(n)%mo1))
                   allocate(agrmet_struc(n)%rlon2_nh4(agrmet_struc(n)%mo1))
                   allocate(agrmet_struc(n)%n112_nh4(agrmet_struc(n)%mo1))
                   call neighbor_interp_input_withgrid(gridDesci4,gridDesco,&
                        agrmet_struc(n)%mo1,agrmet_struc(n)%rlat2_nh4,&
                        agrmet_struc(n)%rlon2_nh4,agrmet_struc(n)%n112_nh4)
!</kluge -- cod testing>
                   elseif(ihemi.eq.2.and.LIS_rc%gridDesc(n,1).eq.0) then
                      allocate(agrmet_struc(n)%rlat1_sh(agrmet_struc(n)%mo2))
                      allocate(agrmet_struc(n)%rlon1_sh(agrmet_struc(n)%mo2))
                      allocate(agrmet_struc(n)%n111_sh(agrmet_struc(n)%mo2))
                      allocate(agrmet_struc(n)%n121_sh(agrmet_struc(n)%mo2))
                      allocate(agrmet_struc(n)%n211_sh(agrmet_struc(n)%mo2))
                      allocate(agrmet_struc(n)%n221_sh(agrmet_struc(n)%mo2))
                      allocate(agrmet_struc(n)%w111_sh(agrmet_struc(n)%mo2))
                      allocate(agrmet_struc(n)%w121_sh(agrmet_struc(n)%mo2))
                      allocate(agrmet_struc(n)%w211_sh(agrmet_struc(n)%mo2))
                      allocate(agrmet_struc(n)%w221_sh(agrmet_struc(n)%mo2))
                      call bilinear_interp_input_withgrid(gridDesci,gridDesco,&
                           agrmet_struc(n)%mo2,agrmet_struc(n)%rlat1_sh,&
                           agrmet_struc(n)%rlon1_sh,agrmet_struc(n)%n111_sh,&
                           agrmet_struc(n)%n121_sh,agrmet_struc(n)%n211_sh,&
                           agrmet_struc(n)%n221_sh,agrmet_struc(n)%w111_sh,&
                           agrmet_struc(n)%w121_sh,agrmet_struc(n)%w211_sh,&
                           agrmet_struc(n)%w221_sh)
                      if(agrmet_struc(n)%imaxnative /= agrmet_struc(n)%imaxgp)then
                         if(agrmet_struc(n)%imax2 == 3600)then
                            gridDesci2(7)=-0.05
                            gridDesci2(3)=900
                         endif
                         allocate(agrmet_struc(n)%rlat1_sh2(agrmet_struc(n)%mo2))
                         allocate(agrmet_struc(n)%rlon1_sh2(agrmet_struc(n)%mo2))
                         allocate(agrmet_struc(n)%n111_sh2(agrmet_struc(n)%mo2))
                         allocate(agrmet_struc(n)%n121_sh2(agrmet_struc(n)%mo2))
                         allocate(agrmet_struc(n)%n211_sh2(agrmet_struc(n)%mo2))
                         allocate(agrmet_struc(n)%n221_sh2(agrmet_struc(n)%mo2))
                         allocate(agrmet_struc(n)%w111_sh2(agrmet_struc(n)%mo2))
                         allocate(agrmet_struc(n)%w121_sh2(agrmet_struc(n)%mo2))
                         allocate(agrmet_struc(n)%w211_sh2(agrmet_struc(n)%mo2))
                         allocate(agrmet_struc(n)%w221_sh2(agrmet_struc(n)%mo2))
                         call bilinear_interp_input_withgrid(gridDesci2,gridDesco,&
                              agrmet_struc(n)%mo2,agrmet_struc(n)%rlat1_sh2,&
                              agrmet_struc(n)%rlon1_sh2,agrmet_struc(n)%n111_sh2,&
                              agrmet_struc(n)%n121_sh2,agrmet_struc(n)%n211_sh2,&
                              agrmet_struc(n)%n221_sh2,agrmet_struc(n)%w111_sh2,&
                              agrmet_struc(n)%w121_sh2,agrmet_struc(n)%w211_sh2,&
                              agrmet_struc(n)%w221_sh2)
                      endif ! nonstandard grid
                      allocate(agrmet_struc(n)%rlat2_sh(agrmet_struc(n)%mo2))
                      allocate(agrmet_struc(n)%rlon2_sh(agrmet_struc(n)%mo2))
                      allocate(agrmet_struc(n)%n112_sh(agrmet_struc(n)%mo2))
                      call neighbor_interp_input_withgrid(gridDesci,gridDesco,&
                           agrmet_struc(n)%mo2,agrmet_struc(n)%rlat2_sh,&
                           agrmet_struc(n)%rlon2_sh,agrmet_struc(n)%n112_sh)
                      if(agrmet_struc(n)%imaxnative /= agrmet_struc(n)%imaxgp)then
                         allocate(agrmet_struc(n)%rlat2_sh2(agrmet_struc(n)%mo2))
                         allocate(agrmet_struc(n)%rlon2_sh2(agrmet_struc(n)%mo2))
                         allocate(agrmet_struc(n)%n112_sh2(agrmet_struc(n)%mo2))
                         call neighbor_interp_input_withgrid(gridDesci2,gridDesco,&
                              agrmet_struc(n)%mo2,agrmet_struc(n)%rlat2_sh2,&
                              agrmet_struc(n)%rlon2_sh2,agrmet_struc(n)%n112_sh2)
                      endif ! nonstandard grid
                      if(agrmet_struc(n)%imaxnative /= agrmet_struc(n)%imaxsmi)then
                         if(agrmet_struc(n)%imax3 == 1440 )then
                            gridDesci3(4)=-0.125
                            gridDesci3(3)=360
                         endif
                         allocate(agrmet_struc(n)%rlat1_sh3(agrmet_struc(n)%mo2))
                         allocate(agrmet_struc(n)%rlon1_sh3(agrmet_struc(n)%mo2))
                         allocate(agrmet_struc(n)%n111_sh3(agrmet_struc(n)%mo2))
                         allocate(agrmet_struc(n)%n121_sh3(agrmet_struc(n)%mo2))
                         allocate(agrmet_struc(n)%n211_sh3(agrmet_struc(n)%mo2))
                         allocate(agrmet_struc(n)%n221_sh3(agrmet_struc(n)%mo2))
                         allocate(agrmet_struc(n)%w111_sh3(agrmet_struc(n)%mo2))
                         allocate(agrmet_struc(n)%w121_sh3(agrmet_struc(n)%mo2))
                         allocate(agrmet_struc(n)%w211_sh3(agrmet_struc(n)%mo2))
                         allocate(agrmet_struc(n)%w221_sh3(agrmet_struc(n)%mo2))
                         call bilinear_interp_input_withgrid(gridDesci3,gridDesco,&
                              agrmet_struc(n)%mo2,agrmet_struc(n)%rlat1_sh3,&
                              agrmet_struc(n)%rlon1_sh3,agrmet_struc(n)%n111_sh3,&
                              agrmet_struc(n)%n121_sh3,agrmet_struc(n)%n211_sh3,&
                              agrmet_struc(n)%n221_sh3,agrmet_struc(n)%w111_sh3,&
                              agrmet_struc(n)%w121_sh3,agrmet_struc(n)%w211_sh3,&
                              agrmet_struc(n)%w221_sh3)
                         allocate(agrmet_struc(n)%rlat2_sh3(agrmet_struc(n)%mo2))
                         allocate(agrmet_struc(n)%rlon2_sh3(agrmet_struc(n)%mo2))
                         allocate(agrmet_struc(n)%n112_sh3(agrmet_struc(n)%mo2))
                         call neighbor_interp_input_withgrid(gridDesci3,gridDesco,&
                              agrmet_struc(n)%mo2,agrmet_struc(n)%rlat2_sh3,&
                              agrmet_struc(n)%rlon2_sh3,agrmet_struc(n)%n112_sh3)
                      endif ! nonstandard grid
!<kluge -- cod testing>
                   allocate(agrmet_struc(n)%rlat1_sh4(agrmet_struc(n)%mo2))
                   allocate(agrmet_struc(n)%rlon1_sh4(agrmet_struc(n)%mo2))
                   allocate(agrmet_struc(n)%n111_sh4(agrmet_struc(n)%mo2))
                   allocate(agrmet_struc(n)%n121_sh4(agrmet_struc(n)%mo2))
                   allocate(agrmet_struc(n)%n211_sh4(agrmet_struc(n)%mo2))
                   allocate(agrmet_struc(n)%n221_sh4(agrmet_struc(n)%mo2))
                   allocate(agrmet_struc(n)%w111_sh4(agrmet_struc(n)%mo2))
                   allocate(agrmet_struc(n)%w121_sh4(agrmet_struc(n)%mo2))
                   allocate(agrmet_struc(n)%w211_sh4(agrmet_struc(n)%mo2))
                   allocate(agrmet_struc(n)%w221_sh4(agrmet_struc(n)%mo2))
                   call bilinear_interp_input_withgrid(gridDesci4,gridDesco,&
                        agrmet_struc(n)%mo2,agrmet_struc(n)%rlat1_sh4,&
                        agrmet_struc(n)%rlon1_sh4,agrmet_struc(n)%n111_sh4,&
                        agrmet_struc(n)%n121_sh4,agrmet_struc(n)%n211_sh4,&
                        agrmet_struc(n)%n221_sh4,agrmet_struc(n)%w111_sh4,&
                        agrmet_struc(n)%w121_sh4,agrmet_struc(n)%w211_sh4,&
                        agrmet_struc(n)%w221_sh4)
                   allocate(agrmet_struc(n)%rlat2_sh4(agrmet_struc(n)%mo2))
                   allocate(agrmet_struc(n)%rlon2_sh4(agrmet_struc(n)%mo2))
                   allocate(agrmet_struc(n)%n112_sh4(agrmet_struc(n)%mo2))
                   call neighbor_interp_input_withgrid(gridDesci4,gridDesco,&
                        agrmet_struc(n)%mo2,agrmet_struc(n)%rlat2_sh4,&
                        agrmet_struc(n)%rlon2_sh4,agrmet_struc(n)%n112_sh4)
!</kluge -- cod testing>
                   endif ! ihemi
                endif ! gridspan 
             enddo

! The below is for non-latlon LIS domain, etc, which has been modified to be
! I left this here for posterity and historical/comparative purposes; obviously the conditional won't be satisfied. - Michael Shaw

          if(LIS_rc%gridDesc(n,1).eq.323232) then
          !elseif(LIS_rc%gridDesc(n,1).ne.0) then

             if(agrmet_struc(n)%imaxnative /= agrmet_struc(n)%imaxgp .or. &
                agrmet_struc(n)%imaxnative /= agrmet_struc(n)%imaxsmi .or. &
                agrmet_struc(n)%imaxnative /= 512)then
                write(*,*) 'Currently mismatched and non 8th mesh '
                write(*,*) 'precip not supported for interpolation'
                write(*,*) 'to non lat-lon LIS grid - stay tuned'
                write(*,*) 'as it is almost there'
                call LIS_endrun()
             endif
 
             allocate(agrmet_struc(n)%smask1(LIS_rc%lnc(n),LIS_rc%lnr(n)))
             allocate(agrmet_struc(n)%smask2(LIS_rc%lnc(n),LIS_rc%lnr(n)))
             agrmet_struc(n)%fillflag1 = .true. 
             agrmet_struc(n)%fillflag2 = .true.
             if(LIS_rc%gridDesc(n,4).ge.0.and.LIS_rc%gridDesc(n,10).ge.0) then 
                agrmet_struc(n)%gridspan = 1 ! NH only
                agrmet_struc(n)%shemi = 1
                agrmet_struc(n)%nhemi = 1
             elseif(LIS_rc%gridDesc(n,4).le.0.and.LIS_rc%gridDesc(n,10).ge.0) then 
                agrmet_struc(n)%gridspan = 2 ! SH only 
                agrmet_struc(n)%shemi = 2
                agrmet_struc(n)%nhemi = 2
             elseif(LIS_rc%gridDesc(n,10).eq.-100.and.LIS_rc%gridDesc(n,20).eq.-100)then 

!-----------------------------------------------------------------------------
! Global grid in polar stereographic projection. No interpolation 
! will be done.
!-----------------------------------------------------------------------------

                agrmet_struc(n)%interp = 0 
                agrmet_struc(n)%gridspan = 3
                agrmet_struc(n)%shemi = 1 
                agrmet_struc(n)%nhemi = 2
             else
                write(*,*) 'Currently spanning across hemispheres is'
                write(*,*) 'not supported'
                call LIS_endrun()
             endif
             agrmet_struc(n)%hemi_nc = LIS_rc%gridDesc(n,2)
             if(agrmet_struc(n)%gridspan.eq.1) then 
                agrmet_struc(n)%hemi_nr(1) = LIS_rc%gridDesc(n,3)
                agrmet_struc(n)%hemi_nr(2) = 0 
                agrmet_struc(n)%mo1 = agrmet_struc(n)%hemi_nc(1)*&
                     agrmet_struc(n)%hemi_nr(1)
                agrmet_struc(n)%mo2 = 0 
             elseif(agrmet_struc(n)%gridspan.eq.2) then 
                agrmet_struc(n)%hemi_nr(1) = 0 
                agrmet_struc(n)%hemi_nr(2) = LIS_rc%gridDesc(n,3)
                agrmet_struc(n)%mo1 = 0 
                agrmet_struc(n)%mo2 = agrmet_struc(n)%hemi_nc(2)*&
                     agrmet_struc(n)%hemi_nr(2)
             else
                agrmet_struc(n)%hemi_nr(1) = LIS_rc%gridDesc(n,3)
                agrmet_struc(n)%hemi_nr(2) = LIS_rc%gridDesc(n,3)
                agrmet_struc(n)%mo1 = agrmet_struc(n)%hemi_nc(1)*&
                     agrmet_struc(n)%hemi_nr(1)
                agrmet_struc(n)%mo2 = agrmet_struc(n)%hemi_nc(2)*&
                     agrmet_struc(n)%hemi_nr(2)             
             endif
             gridDesco = 0 
             gridDesci = 0 
             do ihemi = agrmet_struc(n)%shemi,agrmet_struc(n)%nhemi
                gridDesco(1) = LIS_rc%gridDesc(n,1)
                gridDesco(2) = agrmet_struc(n)%hemi_nc(ihemi)
                gridDesco(3) = agrmet_struc(n)%hemi_nr(ihemi)
                gridDesco(5) = LIS_rc%gridDesc(n,5)
                gridDesco(8) = LIS_rc%gridDesc(n,8)
                gridDesco(6) = LIS_rc%gridDesc(n,6)
                gridDesco(9) = LIS_rc%gridDesc(n,9)
                gridDesco(10) = LIS_rc%gridDesc(n,10)
                gridDesco(11) = LIS_rc%gridDesc(n,11)
                gridDesco(20) = 255
                if(agrmet_struc(n)%gridspan.eq.1.or.&
                     agrmet_struc(n)%gridspan.eq.2) then 
                   gridDesco(4) = LIS_rc%gridDesc(n,4)
                   gridDesco(7) = LIS_rc%gridDesc(n,7)
                endif
                
                if(ihemi.eq.1) then 
                   xmesh = xmeshl
                   orient = 100.0
                else
                   xmesh = -xmeshl
                   orient = 280.0
                endif
                xj1 = float(1)-ypnmcaf
                xi1 = float(1)-xpnmcaf
                
                call polarToLatLon(xi1,xj1,xmesh,orient,alat1,alon1)
                
                gridDesci = 0 
                gridDesci(1) = 5
                gridDesci(2) = agrmet_struc(n)%imax
                gridDesci(3) = agrmet_struc(n)%jmax
                gridDesci(4) = alat1
                gridDesci(5) = alon1
                gridDesci(6) = 8
                gridDesci(7) = orient
                gridDesci(8) = xmesh
                gridDesci(9) = xmesh
                gridDesci(10) = 0.0
                !if(ihemi .eq.2) then 
                   gridDesci(20) = 128
                   gridDesci(11) = 128
                !endif
                
                !gridDesci(20) = 0 
                gridDesci(13) = 1  !global grid
                
                if(agrmet_struc(n)%gridspan.eq.1) then 
                   allocate(agrmet_struc(n)%rlat1_nh(agrmet_struc(n)%mo1))
                   allocate(agrmet_struc(n)%rlon1_nh(agrmet_struc(n)%mo1))
                   allocate(agrmet_struc(n)%n111_nh(agrmet_struc(n)%mo1))
                   allocate(agrmet_struc(n)%n121_nh(agrmet_struc(n)%mo1))
                   allocate(agrmet_struc(n)%n211_nh(agrmet_struc(n)%mo1))
                   allocate(agrmet_struc(n)%n221_nh(agrmet_struc(n)%mo1))
                   allocate(agrmet_struc(n)%w111_nh(agrmet_struc(n)%mo1))
                   allocate(agrmet_struc(n)%w121_nh(agrmet_struc(n)%mo1))
                   allocate(agrmet_struc(n)%w211_nh(agrmet_struc(n)%mo1))
                   allocate(agrmet_struc(n)%w221_nh(agrmet_struc(n)%mo1))
                   call bilinear_interp_input_withgrid(gridDesci,gridDesco,&
                        agrmet_struc(n)%mo1,agrmet_struc(n)%rlat1_nh,&
                        agrmet_struc(n)%rlon1_nh,agrmet_struc(n)%n111_nh,&
                        agrmet_struc(n)%n121_nh,agrmet_struc(n)%n211_nh,&
                        agrmet_struc(n)%n221_nh,agrmet_struc(n)%w111_nh,&
                        agrmet_struc(n)%w121_nh,agrmet_struc(n)%w211_nh,&
                        agrmet_struc(n)%w221_nh)
                   allocate(agrmet_struc(n)%rlat2_nh(agrmet_struc(n)%mo1))
                   allocate(agrmet_struc(n)%rlon2_nh(agrmet_struc(n)%mo1))
                   allocate(agrmet_struc(n)%n112_nh(agrmet_struc(n)%mo1))
                   call neighbor_interp_input_withgrid(gridDesci,gridDesco,&
                        agrmet_struc(n)%mo1,agrmet_struc(n)%rlat2_nh,&
                        agrmet_struc(n)%rlon2_nh,agrmet_struc(n)%n112_nh)
                elseif(agrmet_struc(n)%gridspan.eq.2) then 
                   
                   allocate(agrmet_struc(n)%rlat1_sh(agrmet_struc(n)%mo2))
                   allocate(agrmet_struc(n)%rlon1_sh(agrmet_struc(n)%mo2))
                   allocate(agrmet_struc(n)%n111_sh(agrmet_struc(n)%mo2))
                   allocate(agrmet_struc(n)%n121_sh(agrmet_struc(n)%mo2))
                   allocate(agrmet_struc(n)%n211_sh(agrmet_struc(n)%mo2))
                   allocate(agrmet_struc(n)%n221_sh(agrmet_struc(n)%mo2))
                   allocate(agrmet_struc(n)%w111_sh(agrmet_struc(n)%mo2))
                   allocate(agrmet_struc(n)%w121_sh(agrmet_struc(n)%mo2))
                   allocate(agrmet_struc(n)%w211_sh(agrmet_struc(n)%mo2))
                   allocate(agrmet_struc(n)%w221_sh(agrmet_struc(n)%mo2))
                   call bilinear_interp_input_withgrid(gridDesci,gridDesco,&
                        agrmet_struc(n)%mo2,agrmet_struc(n)%rlat1_sh,&
                        agrmet_struc(n)%rlon1_sh,&
                        agrmet_struc(n)%n111_sh,agrmet_struc(n)%n121_sh,&
                        agrmet_struc(n)%n211_sh,agrmet_struc(n)%n221_sh,&
                        agrmet_struc(n)%w111_sh,agrmet_struc(n)%w121_sh,&
                        agrmet_struc(n)%w211_sh,agrmet_struc(n)%w221_sh)
                   allocate(agrmet_struc(n)%rlat2_sh(agrmet_struc(n)%mo2))
                   allocate(agrmet_struc(n)%rlon2_sh(agrmet_struc(n)%mo2))
                   allocate(agrmet_struc(n)%n112_sh(agrmet_struc(n)%mo2))
                   call neighbor_interp_input_withgrid(gridDesci,gridDesco,&
                        agrmet_struc(n)%mo2,agrmet_struc(n)%rlat2_sh,&
                        agrmet_struc(n)%rlon2_sh,agrmet_struc(n)%n112_sh)
                else

!-------------------------------------------------------------------------
!   No interpolation is being done. So no weights will be computed. 
!-------------------------------------------------------------------------

                   if(ihemi.eq.1) then 
                      allocate(agrmet_struc(n)%rlat1_nh(agrmet_struc(n)%mo1))
                      allocate(agrmet_struc(n)%rlon1_nh(agrmet_struc(n)%mo1))
                      allocate(agrmet_struc(n)%n111_nh(agrmet_struc(n)%mo1))
                      allocate(agrmet_struc(n)%n121_nh(agrmet_struc(n)%mo1))
                      allocate(agrmet_struc(n)%n211_nh(agrmet_struc(n)%mo1))
                      allocate(agrmet_struc(n)%n221_nh(agrmet_struc(n)%mo1))
                      allocate(agrmet_struc(n)%w111_nh(agrmet_struc(n)%mo1))
                      allocate(agrmet_struc(n)%w121_nh(agrmet_struc(n)%mo1))
                      allocate(agrmet_struc(n)%w211_nh(agrmet_struc(n)%mo1))
                      allocate(agrmet_struc(n)%w221_nh(agrmet_struc(n)%mo1))
                      allocate(agrmet_struc(n)%rlat2_nh(agrmet_struc(n)%mo1))
                      allocate(agrmet_struc(n)%rlon2_nh(agrmet_struc(n)%mo1))
                      allocate(agrmet_struc(n)%n112_nh(agrmet_struc(n)%mo1))
                   elseif(ihemi.eq.2) then 
                      allocate(agrmet_struc(n)%rlat1_sh(agrmet_struc(n)%mo2))
                      allocate(agrmet_struc(n)%rlon1_sh(agrmet_struc(n)%mo2))
                      allocate(agrmet_struc(n)%n111_sh(agrmet_struc(n)%mo2))
                      allocate(agrmet_struc(n)%n121_sh(agrmet_struc(n)%mo2))
                      allocate(agrmet_struc(n)%n211_sh(agrmet_struc(n)%mo2))
                      allocate(agrmet_struc(n)%n221_sh(agrmet_struc(n)%mo2))
                      allocate(agrmet_struc(n)%w111_sh(agrmet_struc(n)%mo2))
                      allocate(agrmet_struc(n)%w121_sh(agrmet_struc(n)%mo2))
                      allocate(agrmet_struc(n)%w211_sh(agrmet_struc(n)%mo2))
                      allocate(agrmet_struc(n)%w221_sh(agrmet_struc(n)%mo2))
                      allocate(agrmet_struc(n)%rlat2_sh(agrmet_struc(n)%mo2))
                      allocate(agrmet_struc(n)%rlon2_sh(agrmet_struc(n)%mo2))
                      allocate(agrmet_struc(n)%n112_sh(agrmet_struc(n)%mo2))
                   endif
                endif
             enddo
          endif
          
! Read the precip climatologies and set alarms

          agrmet_struc(n)%pcpclimoAlarmTime = 0.0

          allocate(agrmet_struc(n)%cliprc1(LIS_rc%lnc(n), LIS_rc%lnr(n)))
          agrmet_struc(n)%cliprc1 = 0
          allocate(agrmet_struc(n)%clippd1(LIS_rc%lnc(n), LIS_rc%lnr(n)))
          agrmet_struc(n)%clippd1 = 0
          allocate(agrmet_struc(n)%clirtn1(LIS_rc%lnc(n), LIS_rc%lnr(n)))
          agrmet_struc(n)%clirtn1 = 0
          allocate(agrmet_struc(n)%cliprc2(LIS_rc%lnc(n), LIS_rc%lnr(n)))
          agrmet_struc(n)%cliprc2 = 0
          allocate(agrmet_struc(n)%clippd2(LIS_rc%lnc(n), LIS_rc%lnr(n)))
          agrmet_struc(n)%clippd2 = 0
          allocate(agrmet_struc(n)%clirtn2(LIS_rc%lnc(n), LIS_rc%lnr(n)))
          agrmet_struc(n)%clirtn2 = 0
          
          allocate(agrmet_struc(n)%cliprc(LIS_rc%lnc(n), LIS_rc%lnr(n)))
          agrmet_struc(n)%cliprc = 0
          allocate(agrmet_struc(n)%clippd(LIS_rc%lnc(n), LIS_rc%lnr(n)))
          agrmet_struc(n)%clippd = 0
          allocate(agrmet_struc(n)%clirtn(LIS_rc%lnc(n), LIS_rc%lnr(n)))
          agrmet_struc(n)%clirtn = 0
          
          allocate(agrmet_struc(n)%sfcprs(6,LIS_rc%lnc(n), LIS_rc%lnr(n)))
          agrmet_struc(n)%sfcprs = 0
          allocate(agrmet_struc(n)%sfcrlh(6,LIS_rc%lnc(n), LIS_rc%lnr(n)))
          agrmet_struc(n)%sfcrlh = 0
          allocate(agrmet_struc(n)%sfcspd(6,LIS_rc%lnc(n), LIS_rc%lnr(n)))
          agrmet_struc(n)%sfcspd = 0
          allocate(agrmet_struc(n)%sfctmp(6,LIS_rc%lnc(n), LIS_rc%lnr(n)))
          agrmet_struc(n)%sfctmp = 0
          allocate(agrmet_struc(n)%lasprs(LIS_rc%lnc(n), LIS_rc%lnr(n)))
          agrmet_struc(n)%lasprs = 0
          allocate(agrmet_struc(n)%lasrlh(LIS_rc%lnc(n), LIS_rc%lnr(n)))
          agrmet_struc(n)%lasrlh = 0
          allocate(agrmet_struc(n)%lasspd(LIS_rc%lnc(n), LIS_rc%lnr(n)))
          agrmet_struc(n)%lasspd = 0
          allocate(agrmet_struc(n)%lastmp(LIS_rc%lnc(n), LIS_rc%lnr(n)))
          agrmet_struc(n)%lastmp = 0

          ! EMK BEGIN
!<rm -- jim merge>
#if 0
          allocate(agrmet_struc(n)%sfcprs_glb(6,LIS_rc%gnc(n), LIS_rc%gnr(n)))
#endif
!<rm -- jim merge>
          allocate(agrmet_struc(n)%sfcrlh_glb(6,LIS_rc%gnc(n), LIS_rc%gnr(n)))
          allocate(agrmet_struc(n)%sfcspd_glb(6,LIS_rc%gnc(n), LIS_rc%gnr(n)))
          allocate(agrmet_struc(n)%sfctmp_glb(6,LIS_rc%gnc(n), LIS_rc%gnr(n)))

!<rm -- jim merge>
#if 0
          allocate(agrmet_struc(n)%lasprs_glb(LIS_rc%gnc(n), LIS_rc%gnr(n)))
          allocate(agrmet_struc(n)%lasrlh_glb(LIS_rc%gnc(n), LIS_rc%gnr(n)))
          allocate(agrmet_struc(n)%lasspd_glb(LIS_rc%gnc(n), LIS_rc%gnr(n)))
          allocate(agrmet_struc(n)%lastmp_glb(LIS_rc%gnc(n), LIS_rc%gnr(n)))
#endif
!<rm -- jim merge>
          ! EMK END

! Initialize clippd.  This array is set only when cdfs2swch is 1,
! but it is used when either cdfs2swch or clswch is 1.
! Values greater than -9990 are considered valid.
          agrmet_struc(n)%clippd = -9999.0

! 2 hemispheres, max 4 3 hr intervals. 

          allocate(agrmet_struc(n)%mrgp(LIS_rc%lnc(n), LIS_rc%lnr(n),4))
          agrmet_struc(n)%pcp_start = .true. 
          agrmet_struc(n)%pcp_ready = .false. 
          agrmet_struc(n)%lastSfcalcHour = 0
          agrmet_struc(n)%lastPcpHour = 0
          call AGRMET_read_pcpclimodata(n)

          agrmet_struc(n)%albAlarmTime = 0.0
          agrmet_struc(n)%gfracAlarmTime = 0.0

          if ( ( agrmet_struc(n)%first_guess_source /= 'GFS' ) .and. &
               ( agrmet_struc(n)%first_guess_source /= 'GALWEM' ) ) then
             write(LIS_logunit,*) &
                '[ERR] First guess source is not correctly defined.'
             call LIS_endrun
          endif

          if ( agrmet_struc(n)%galwemprecswch == 1 .and. &
               agrmet_struc(n)%gfsprecswch == 1 ) then
             write(LIS_logunit,*) '[ERR] Cannot enable both ' // &
                                  '"AGRMET use GFS precip:" and ' // &
                                  '"AGRMET use GALWEM precip:".'
             call LIS_endrun
          elseif ( agrmet_struc(n)%first_guess_source == 'GALWEM' .and. &
                   agrmet_struc(n)%gfsprecswch == 1 ) then
             write(LIS_logunit,*) &
                '[ERR] Cannot enable "AGRMET use GFS precip:" ' // &
                'with GALWEM first guess'
             call LIS_endrun
          elseif ( agrmet_struc(n)%first_guess_source == 'GFS' .and. &
                   agrmet_struc(n)%galwemprecswch == 1 ) then
             write(LIS_logunit,*) &
                '[ERR] Cannot enable "AGRMET use GALWEM precip:" ' // &
                'with GFS first guess'
             call LIS_endrun
          endif

          !interpolation weights for conversion from GFS to LIS grid
          gridDesci_glb = 0
          gridDesci_glb(1) = 0
          gridDesci_glb(2) = 720
          gridDesci_glb(3) = 361
          gridDesci_glb(4) = -90.000
          gridDesci_glb(5) = -180.000
          gridDesci_glb(6) = 128
          gridDesci_glb(7) = 90.000
!          gridDesci_glb(8) = 180.000
          gridDesci_glb(8) = 179.500 !EMK TEST BUG FIX...No repeating meridian
          gridDesci_glb(9) = 0.500
          gridDesci_glb(10) = 0.5000
          gridDesci_glb(20) = 0

          call gfs_reset_interp_input(n, findex, gridDesci_glb)

          !interpolation weights for conversion from GALWEM to LIS grid
          ! EMK...Only support 17-km or 10-km resolutions.
          if (agrmet_struc(n)%galwem_res == 17) then
             gridDesci_glb = 0
             gridDesci_glb(1) = 0
             gridDesci_glb(2) = 1536
             gridDesci_glb(3) = 1152
             gridDesci_glb(4) = -89.9219
             gridDesci_glb(5) = -179.882813
             gridDesci_glb(6) = 128
             gridDesci_glb(7) = 89.9219
             gridDesci_glb(8) = 179.887
             gridDesci_glb(9) = 0.234378
             gridDesci_glb(10) = 0.15625
             gridDesci_glb(20) = 0
          else if (agrmet_struc(n)%galwem_res == 10) then
             gridDesci_glb = 0
             gridDesci_glb(1) = 0
             gridDesci_glb(2) = 2560
             gridDesci_glb(3) = 1920
             gridDesci_glb(4) = -89.9531250
             gridDesci_glb(5) = -179.9296875
             gridDesci_glb(6) = 128
             gridDesci_glb(7) = 89.9531250
             gridDesci_glb(8) = 179.9296875
             gridDesci_glb(9) = 0.140625
             gridDesci_glb(10) = 0.093750
             gridDesci_glb(20) = 0
          else
             write(LIS_logunit,*) '[ERR] LIS only supports GALWEM 17 or 10 km'
             call LIS_endrun()
          end if

          call galwem_reset_interp_input(n, findex, gridDesci_glb)

!interpolation weights for conversion from CMORPH to LIS grid

          gridDesciCmor = 0
          gridDesciCmor(1) = 0
          gridDesciCmor(2) = agrmet_struc(n)%imaxcmor
          gridDesciCmor(3) = agrmet_struc(n)%jmaxcmor
          gridDesciCmor(4) = agrmet_struc(n)%cmorminlat
          gridDesciCmor(5) = agrmet_struc(n)%cmorminlon
          gridDesciCmor(6) = 128
          gridDesciCmor(7) = agrmet_struc(n)%cmormaxlat
          gridDesciCmor(8) = agrmet_struc(n)%cmormaxlon
          gridDesciCmor(9) = agrmet_struc(n)%cmordy
          gridDesciCmor(10) = agrmet_struc(n)%cmordx
          gridDesciCmor(20) = 64

          allocate(agrmet_struc(n)%n112cmor(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(agrmet_struc(n)%n122cmor(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(agrmet_struc(n)%n212cmor(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(agrmet_struc(n)%n222cmor(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(agrmet_struc(n)%w112cmor(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(agrmet_struc(n)%w122cmor(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(agrmet_struc(n)%w212cmor(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(agrmet_struc(n)%w222cmor(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
       
          call conserv_interp_input(n,&
               gridDesciCmor,&
               agrmet_struc(n)%n112cmor,agrmet_struc(n)%n122cmor,agrmet_struc(n)%n212cmor,&
               agrmet_struc(n)%n222cmor,agrmet_struc(n)%w112cmor,agrmet_struc(n)%w122cmor,&
               agrmet_struc(n)%w212cmor,agrmet_struc(n)%w222cmor)

       else
          allocate(agrmet_struc(n)%smask1(LIS_rc%lnc(n),LIS_rc%lnr(n)))
          allocate(agrmet_struc(n)%smask2(LIS_rc%lnc(n),LIS_rc%lnr(n)))
          agrmet_struc(n)%fillflag1 = .true. 
          agrmet_struc(n)%fillflag2 = .true.

          !EMK...Original code
!           gridDesci     = 0
!           gridDesci(1)  = 0 
!           gridDesci(2)  = 1440 
!           gridDesci(3)  = 600
!           gridDesci(4)  = -59.875
!           gridDesci(5)  = -179.875
!           gridDesci(6)  = 128
!           gridDesci(7)  = 89.875
!           gridDesci(8)  = 179.875
!           gridDesci(9)  = 0.25
!           gridDesci(10) = 0.25
!           gridDesci(20) = 64

!           agrmet_struc(n)%ncol = 1440 
!           agrmet_struc(n)%nrow = 600

          ! EMK NEW Code...Use new lis.config settings
          gridDesci = 0
          gridDesci(1)  = 0 ! LatLon
          gridDesci(2)  = 1 + &
               int((agrmet_struc(n)%forcingUpperRightLon - &
               agrmet_struc(n)%forcingLowerLeftLon) / &
               agrmet_struc(n)%forcingResDx) ! Columns
          gridDesci(3)  = 1 + &
               int((agrmet_struc(n)%forcingUpperRightLat - &
               agrmet_struc(n)%forcingLowerLeftLat) / &
               agrmet_struc(n)%forcingResDy) ! Rows
          gridDesci(4)  = agrmet_struc(n)%forcingLowerLeftLat
          gridDesci(5)  = agrmet_struc(n)%forcingLowerLeftLon
          gridDesci(6)  = 128   ! Not used
          gridDesci(7)  = agrmet_struc(n)%forcingUpperRightLat
          gridDesci(8)  = agrmet_struc(n)%forcingUpperRightLon
          gridDesci(9)  = agrmet_struc(n)%forcingResDx
          gridDesci(10) = agrmet_struc(n)%forcingResDy
          gridDesci(20) = 64 ! E-W ordering

          agrmet_struc(n)%ncol = gridDesci(2)
          agrmet_struc(n)%nrow = gridDesci(3)
          
          if(trim(LIS_rc%met_interp(findex)).eq."bilinear") then 

             allocate(agrmet_struc(n)%n11_anl(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
             allocate(agrmet_struc(n)%n12_anl(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
             allocate(agrmet_struc(n)%n21_anl(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
             allocate(agrmet_struc(n)%n22_anl(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
             allocate(agrmet_struc(n)%w11_anl(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
             allocate(agrmet_struc(n)%w12_anl(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
             allocate(agrmet_struc(n)%w21_anl(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
             allocate(agrmet_struc(n)%w22_anl(LIS_rc%lnc(n)*LIS_rc%lnr(n)))

             call bilinear_interp_input(n,gridDesci(:),&
                  agrmet_struc(n)%n11_anl,agrmet_struc(n)%n12_anl,&
                  agrmet_struc(n)%n21_anl,agrmet_struc(n)%n22_anl,&
                  agrmet_struc(n)%w11_anl,agrmet_struc(n)%w12_anl,&
                  agrmet_struc(n)%w21_anl,agrmet_struc(n)%w22_anl)
          else if(trim(LIS_rc%met_interp(findex)).eq."average") then
             agrmet_struc(n)%mi111 = agrmet_struc(n)%nrow * &
                  agrmet_struc(n)%ncol
             allocate(agrmet_struc(n)%n111_anl(agrmet_struc(n)%mi111))
             call upscaleByAveraging_input(gridDesci, &
                  LIS_rc%gridDesc(n,:), agrmet_struc(n)%mi111, &
                  LIS_rc%lnc(n)*LIS_rc%lnr(n), &
                  agrmet_struc(n)%n111_anl)
          else
             write(LIS_logunit,*) &
                  '[ERR] Unsupported AGRMET interpolation option ', &
                  trim(LIS_rc%met_interp(findex))
             write(LIS_logunit,*) '[ERR] Aborting...'
             call LIS_endrun()
          endif
       endif

       if(LIS_rc%pcp_downscale(findex).ne.0) then
          call LIS_init_pcpclimo_native(n,findex,&
               agrmet_struc(n)%ncol, &  
               agrmet_struc(n)%nrow )
          
       endif

    enddo

    call AGRMET_sanity_check

  end subroutine init_AGRMET

  subroutine AGRMET_sanity_check()

     use LIS_coreMod,    only : LIS_rc
     use LIS_albedoMod,  only : LIS_alb
     use LIS_vegDataMod, only : LIS_gfrac
     use LIS_snowMod,    only : LIS_snow_struc
     use LIS_logMod

     implicit none

     integer :: n

     do n = 1, LIS_rc%nnest
        if ( .not. allocated(LIS_alb(n)%albsf) ) then
           write(LIS_logunit,*) &
              "[ERR] Required variable LIS_alb%albsf is not allocated."
           call LIS_endrun
        endif

        if ( .not. allocated(LIS_gfrac(n)%greenness) ) then
           write(LIS_logunit,*) &
              "[ERR] Required variable LIS_gfrac%greenness is not allocated."
           call LIS_endrun
        endif

        if ( .not. allocated(LIS_alb(n)%mxsnalb) ) then
           write(LIS_logunit,*) &
              "[ERR] Required variable LIS_alb%mxsnalb is not allocated."
           call LIS_endrun
        endif

        if ( .not. allocated(LIS_snow_struc(n)%sneqv) ) then
           write(LIS_logunit,*) &
              "[ERR] Required variable LIS_snow_struc%sneqv is not allocated."
           call LIS_endrun
        endif
     enddo
  end subroutine AGRMET_sanity_check

end module AGRMET_forcingMod


