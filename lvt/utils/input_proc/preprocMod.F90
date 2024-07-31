!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
module preprocMod
  use ESMF
  use netcdf
  use map_utils

  public :: LVT_rc
  public :: LVT_config
  public :: LVT_domain
  public :: LVT_LSMparam_struc  

  public :: LVT_readData   
  public :: ParamProcInit
  public :: paramProcWrite
  public :: LVT_verify

  type lvtrcdec

! -- Tile parameters:
     character*50           :: lis_map_proj
     character*50           :: vegsrc
     integer                :: nt          ! Number veg types 
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
     
     real                   :: udef

     integer                :: nnest
     integer, allocatable       :: gnc(:)
     integer, allocatable       :: gnr(:)
     integer, allocatable       :: lnc(:)
     integer, allocatable       :: lnr(:)
     integer, allocatable       :: lnc_red(:)
     integer, allocatable       :: lnr_red(:)


! -- Parameter grid description arrays:
     real, allocatable          :: gridDesc(:,:)
     real, allocatable          :: lc_gridDesc(:,:)
     character*50           :: lc_proj
     character*50,  allocatable :: mask_type(:)

! -- Parameter filepath names:
     character*140, allocatable :: mfile(:)
     character*140, allocatable :: vfile(:)
     integer,       allocatable :: vfile_form(:)

     integer                :: bareclass 
     integer                :: urbanclass
     integer                :: snowclass 
     integer                :: waterclass
     integer                :: lakeclass
     integer                :: wetlandclass
     integer                :: glacierclass
     integer                :: permafrostclass
  
  end type lvtrcdec

  type, public :: lvt_domain_type 
     type(proj_info)        :: lvtproj
     type(proj_info)        :: lisproj
     integer, allocatable       :: gindex(:,:)
     integer, allocatable       :: lis_gindex(:,:)
     integer, allocatable       :: ntiles_pergrid(:)
     integer, allocatable       :: str_tind(:) !starting tile id for each grid
     real                   :: minLat, maxLat, minLon, maxLon
  end type lvt_domain_type
  
  type, public :: LVT_paramEntry
     character*50  :: short_name
     integer       :: selectOpt
     character*50  :: source
     character*20  :: units
     integer       :: vid
     integer       :: vlevels
     integer       :: num_times
     integer       :: num_bins
     real          :: valid_min
     real          :: valid_max
     character*100 :: standard_name
     real, allocatable :: value(:,:,:)
  end type LVT_paramEntry

!- LSM-specific parameters:
  type, public :: lsmparam_type_dec
     character*100 :: param_filename
     
     type(LVT_paramEntry) :: landmask
     type(LVT_paramEntry) :: landcover
     type(LVT_paramEntry) :: sfctype     ! Surface model type

     type(LVT_paramEntry) :: texture     ! Soil texture map
     type(LVT_paramEntry) :: sand
     type(LVT_paramEntry) :: clay
     type(LVT_paramEntry) :: silt
     type(LVT_paramEntry) :: soilsfgrd   ! Soils tile gridcell fraction

     type(LVT_paramEntry) :: elevation
     type(LVT_paramEntry) :: elevfgrd    ! Elev tile gridcell fraction 
     type(LVT_paramEntry) :: slope
     type(LVT_paramEntry) :: slopefgrd   ! Slope tile gridcell fraction
     type(LVT_paramEntry) :: aspect
     type(LVT_paramEntry) :: aspectfgrd  ! Aspect tile gridcell fraction
     
  end type lsmparam_type_dec


  type(lvtrcdec)    :: LVT_rc
  type(ESMF_Config) :: LVT_config
  type(lvt_domain_type),allocatable :: LVT_domain(:)
  type(lsmparam_type_dec), allocatable :: LVT_LSMparam_struc(:)

!BOP
! 
!  !ROUTINE: LVT_readData
! \label{LVT_readData}
! 
! !INTERFACE: 
  interface LVT_readData

! !PRIVATE MEMBER FUNCTIONS: 
     module procedure read2Ddata
     module procedure readLCdata
! 
! !DESCRIPTION: 
!  Routine to read 2d data from a binary file, with direct access format. A special
!  routine is required to read the landcover data since it includes a distribution of 
!  vegetation types at each grid point. 
!  
!EOP
  end interface

contains

!BOP
! !ROUTINE: paramProcInit
! \label{paramProcInit}
!
! !INTERFACE: 
  subroutine paramProcInit()

! !USES:

    integer   :: n 
    integer   :: rc
    integer   :: k

    allocate(LVT_LSMparam_struc(LVT_rc%nnest))

    call ESMF_ConfigFindLabel(LVT_config,&
         "Processed LSM parameter filename:",rc=rc)
    call ESMF_ConfigGetAttribute(LVT_config,&
         LVT_LSMparam_struc(1)%param_filename,rc=rc)
    call LVT_verify(rc,"Processed LSM parameter filename: not defined")

    call ESMF_ConfigFindLabel(LVT_config,"Landcover data source:",rc=rc)
    call ESMF_ConfigGetAttribute(LVT_config,LVT_rc%vegsrc,rc=rc)
    call LVT_verify(rc,"Landcover data source: not defined")

! - LSM-based parameters:

  !- Set land cover, mask and surface type sources:
     LVT_LSMparam_struc(1)%landcover%source  = trim(LVT_rc%vegsrc)
     LVT_LSMparam_struc(1)%sfctype%source    = trim(LVT_rc%vegsrc)
     LVT_LSMparam_struc(1)%landmask%source   = trim(LVT_rc%vegsrc)
     LVT_LSMparam_struc(1)%landcover%standard_name= trim(LVT_rc%vegsrc)//" land cover"
     LVT_LSMparam_struc(1)%landmask%standard_name = trim(LVT_rc%vegsrc)//" land mask"
     LVT_LSMparam_struc(1)%sfctype%standard_name  = trim(LVT_rc%vegsrc)//" surface type"

  !- Land cover map:
     if( LVT_rc%vegsrc == "UMD" ) then
       LVT_LSMparam_struc(1)%landcover%num_bins   = 13
       LVT_LSMparam_struc(1)%landcover%vlevels    = 13
       LVT_LSMparam_struc(1)%sfctype%num_bins     = 13
       LVT_LSMparam_struc(1)%sfctype%vlevels      = 13

     elseif( LVT_rc%vegsrc == "USGS" ) then
       LVT_LSMparam_struc(1)%landcover%num_bins   = 24 
       LVT_LSMparam_struc(1)%landcover%vlevels    = 24
       LVT_LSMparam_struc(1)%sfctype%num_bins     = 24
       LVT_LSMparam_struc(1)%sfctype%vlevels      = 24

     elseif( LVT_rc%vegsrc == "MODIS" ) then
       LVT_LSMparam_struc(1)%landcover%num_bins   = 20
       LVT_LSMparam_struc(1)%landcover%vlevels    = 20
       LVT_LSMparam_struc(1)%sfctype%num_bins     = 20
       LVT_LSMparam_struc(1)%sfctype%vlevels      = 20
     end if

     LVT_LSMparam_struc(1)%landcover%short_name = "LANDCOVER"
     LVT_LSMparam_struc(1)%landcover%selectOpt  = 1
     LVT_LSMparam_struc(1)%landcover%units      = "-"
     LVT_LSMparam_struc(1)%landcover%num_times  = 1

     LVT_LSMparam_struc(1)%sfctype%short_name = "SURFACETYPE"
     LVT_LSMparam_struc(1)%sfctype%selectOpt = 1
     LVT_LSMparam_struc(1)%sfctype%units     = "-"
     LVT_LSMparam_struc(1)%sfctype%num_times = 1

     LVT_LSMparam_struc(1)%landmask%short_name = "LANDMASK"
     LVT_LSMparam_struc(1)%landmask%num_bins   = 1
     LVT_LSMparam_struc(1)%landmask%vlevels    = 1
     LVT_LSMparam_struc(1)%landmask%selectOpt  = 1
     LVT_LSMparam_struc(1)%landmask%units      = "-"
     LVT_LSMparam_struc(1)%landmask%num_times  = 1

 !- Set number parameter layers based on "vertical" levels (set in param_attribs.txt):
    LVT_rc%nt = LVT_LSMparam_struc(1)%landcover%num_bins

!-- Initialize and/or read in each parameter selected in preproc.config:
    call LMLC_init

  end subroutine ParamProcInit

!BOP
! 
! !ROUTINE: LMLC_init
! \label{LMLC_init}
! 
! !INTERFACE:
  subroutine LMLC_init()

! !USES:
    !NONE
! 
! !DESCRIPTION:
!
! Allocates memory for data structures for reading 
! landmask and landcover datasets
!
!EOP
    implicit none
    integer :: rc
    integer :: n
! _________________________________________________

  ! Initialize gridcell water fraction (default value 0.5)
    allocate(LVT_rc%lc_gridDesc(LVT_rc%nnest,20))
    
!-- Allocate landcover and landmask file arrays:
    do n = 1, LVT_rc%nnest
       allocate(LVT_LSMparam_struc(n)%landmask%value(&
                LVT_rc%lnc(n),LVT_rc%lnr(n),&
               (LVT_LSMparam_struc(n)%landmask%num_bins)))

       allocate(LVT_LSMparam_struc(n)%landcover%value(&
                LVT_rc%lnc(n),LVT_rc%lnr(n),&
                LVT_LSMparam_struc(n)%landcover%num_bins))

       allocate(LVT_LSMparam_struc(n)%sfctype%value(&
                LVT_rc%lnc(n),LVT_rc%lnr(n),&
                LVT_LSMparam_struc(n)%sfctype%num_bins))

   !- Read in options from input config file:
      call ESMF_ConfigFindLabel(LVT_config,"landmask file:",rc=rc)
      call ESMF_ConfigGetAttribute(LVT_config,LVT_rc%mfile(n),rc=rc)
      call LVT_verify(rc,'landmask file: not specified')

!      if( LVT_rc%gridDesc(n,9) .ne. LVT_rc%mask_gridDesc(n,9) ) then 
!         write(*,*) "For now, the landmask resolution must be same as "
!         write(*,*) " the LIS-run domain.  Please either provide a landmask"
!         write(*,*) " that matches the run-domain or select 'create' mask option."
!         write(*,*) " Stopping ..."
!         call LVT_endrun
!      endif
    enddo

!-- Landcover dataset file and option inputs:

    call ESMF_ConfigFindLabel(LVT_config,"landcover file:",rc=rc)
    do n=1,LVT_rc%nnest
       call ESMF_ConfigGetAttribute(LVT_config,LVT_rc%vfile(n),rc=rc)
       call LVT_verify(rc,'landcover file: not specified')
    enddo

    call ESMF_ConfigFindLabel(LVT_config,"landcover file format:",rc=rc)
    do n=1,LVT_rc%nnest
       call ESMF_ConfigGetAttribute(LVT_config,LVT_rc%vfile_form(n),rc=rc)
       call LVT_verify(rc,'landcover file format: not specified')
    enddo

    call ESMF_ConfigGetAttribute(LVT_config,LVT_rc%lc_proj,&
         label="landcover map projection:",rc=rc)
    call LVT_verify(rc,'landcover map projection: option not specified in the config file')

    call readDomainConfigSpecs("landcover", LVT_rc%lc_proj, LVT_rc%lc_gridDesc)

       
!-- Read land cover and mask files:
    do n=1,LVT_rc%nnest

!       write(*,*) "Reading landmask values"
       if( LVT_rc%vegsrc == "UMD" ) then
          call read_maskfile( n, LVT_LSMparam_struc(n)%landmask%value )
          
          !       write(*,*) "Reading landcover values"
          call read_UMDlc(n, &
               LVT_LSMparam_struc(n)%landcover%num_bins, &
               LVT_LSMparam_struc(n)%landcover%value,       &
               LVT_LSMparam_struc(n)%landmask%value,        &
               LVT_LSMparam_struc(n)%sfctype%value  )
       elseif(LVT_rc%vegsrc == "USGS") then 
          call read_maskfile( n, LVT_LSMparam_struc(n)%landmask%value )
          
          !       write(*,*) "Reading landcover values"
          call read_USGSlc(n, &
               LVT_LSMparam_struc(n)%landcover%num_bins, &
               LVT_LSMparam_struc(n)%landcover%value,       &
               LVT_LSMparam_struc(n)%landmask%value,        &
               LVT_LSMparam_struc(n)%sfctype%value  )

       elseif(LVT_rc%vegsrc == "MODIS") then 
          call read_maskfile( n, LVT_LSMparam_struc(n)%landmask%value )
          
          !       write(*,*) "Reading landcover values"
          call read_IGBP_MODISlc(n, &
               LVT_LSMparam_struc(n)%landcover%num_bins, &
               LVT_LSMparam_struc(n)%landcover%value,       &
               LVT_LSMparam_struc(n)%landmask%value,        &
               LVT_LSMparam_struc(n)%sfctype%value  )
          
       endif
       
          
       LVT_LSMparam_struc(n)%landcover%vlevels = &
            LVT_LSMparam_struc(n)%landcover%num_bins
       LVT_LSMparam_struc(n)%sfctype%vlevels = &
            LVT_LSMparam_struc(n)%sfctype%num_bins

    enddo

    write(*,*) "Finished reading landmask and landcover data"

  end subroutine LMLC_init

subroutine read_maskfile(nest, localmask)

! !USES:

  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: nest
  real                :: localmask(LVT_rc%lnc(nest),LVT_rc%lnr(nest))

! !DESCRIPTION:
!  This subroutine reads the landmask data and returns the 
!  array in a lat/lon projection.   
!
!  The arguments are:
!  \begin{description}
!   \item[n]
!    index of nest
!   \item[localmask]
!    landmask for the region of interest
!   \end{description}
!
!EOP      

  integer :: ios1
  integer :: ftn
  logical :: file_exists

  inquire(file=trim(LVT_rc%mfile(nest)), exist=file_exists)

  if(file_exists) then 
     write(*,*)'MSG: LAT/LON mask -- Reading ',trim(LVT_rc%mfile(nest))
     ftn = 100
     open(ftn,file=trim(LVT_rc%mfile(nest)),form='unformatted',recl=4, & 
          access='direct',iostat=ios1)
     
     call LVT_readData(nest, ftn, LVT_rc%lc_proj, LVT_rc%lc_gridDesc(nest,:), localmask)

#if ( defined INC_WATER_PTS )
     localmask = 1.0
#endif

     close(ftn)

  else
     write(*,*) 'landmask map: ',trim(LVT_rc%mfile(nest)), ' does not exist'
     write(*,*) 'program stopping ...'
     stop
  endif

end subroutine read_maskfile

!BOP
! !ROUTINE: readDomainConfigSpecs
! \label{readDomainConfigSpecs}
! 
! !INTERFACE: 
 subroutine readDomainConfigSpecs(segment_name, param_proj, domain_info)

! !USES: 
   !NONE

! !ARGUMENTS: 
   character(len=*),  intent(in)    :: segment_name
   character(50),     intent(in)    :: param_proj
   real,              intent(inout) :: domain_info(LVT_rc%nnest,20)
! 
! !DESCRIPTION: 
!   This subroutine reads the relevant section of the lis.config file 
!   that defines the domain specifications for a particular surface 
!   parameter dataset, based on the map projection used to define 
!   the surface parameter dataset. 
!EOP
   integer  :: i, rc

! ---
  select case ( param_proj )

    case ( "latlon" )  

      call ESMF_ConfigFindLabel(LVT_config,trim(segment_name)//" lower left lat:",rc=rc)
      do i=1,LVT_rc%nnest
         call ESMF_ConfigGetAttribute(LVT_config,domain_info(i,4),rc=rc)
         call LVT_verify(rc,'please specify '//trim(segment_name)//' lower left lat:')
      enddo
      
      call ESMF_ConfigFindLabel(LVT_config,trim(segment_name)//" lower left lon:",rc=rc)
      do i=1,LVT_rc%nnest
         call ESMF_ConfigGetAttribute(LVT_config,domain_info(i,5),rc=rc)
         call LVT_verify(rc,'please specify '//trim(segment_name)//' lower left lon:')
      enddo
      
      call ESMF_ConfigFindLabel(LVT_config,trim(segment_name)//" upper right lat:",rc=rc)
      do i=1,LVT_rc%nnest
         call ESMF_ConfigGetAttribute(LVT_config,domain_info(i,7),rc=rc)
         call LVT_verify(rc,'please specify '//trim(segment_name)//' upper right lat:')
      enddo
      
      call ESMF_ConfigFindLabel(LVT_config,trim(segment_name)//" upper right lon:",rc=rc)
      do i=1,LVT_rc%nnest
         call ESMF_ConfigGetAttribute(LVT_config,domain_info(i,8),rc=rc)
         call LVT_verify(rc,'please specify '//trim(segment_name)//' upper right lon:')
      enddo
      
      call ESMF_ConfigFindLabel(LVT_config,trim(segment_name)//" resolution (dx):",rc=rc)
      do i=1,LVT_rc%nnest
         call ESMF_ConfigGetAttribute(LVT_config,domain_info(i,9),rc=rc)
         call LVT_verify(rc,'please specify '//trim(segment_name)//' resolution (dx):')
      enddo
      
      call ESMF_ConfigFindLabel(LVT_config,trim(segment_name)//" resolution (dy):",rc=rc)
      do i=1,LVT_rc%nnest
         call ESMF_ConfigGetAttribute(LVT_config,domain_info(i,10),rc=rc)
         call LVT_verify(rc,'please specify '//trim(segment_name)//' resolution (dy):')
      enddo
      
! ---
    case ( "gaussian" )  

      call ESMF_ConfigFindLabel(LVT_config,trim(segment_name)//" first grid point lat:",rc=rc)
      do i=1,LVT_rc%nnest
         call ESMF_ConfigGetAttribute(LVT_config,domain_info(i,1),rc=rc)
         call LVT_verify(rc,'please specify '//trim(segment_name)//' first grid point lat:')
      enddo
      
      call ESMF_ConfigFindLabel(LVT_config,trim(segment_name)//" first grid point lon:",rc=rc)
      do i=1,LVT_rc%nnest
         call ESMF_ConfigGetAttribute(LVT_config,domain_info(i,2),rc=rc)
         call LVT_verify(rc,'please specify '//trim(segment_name)//' first grid point lon:')
      enddo
      
      call ESMF_ConfigFindLabel(LVT_config,trim(segment_name)//" last grid point lat:",rc=rc)
      do i=1,LVT_rc%nnest
         call ESMF_ConfigGetAttribute(LVT_config,domain_info(i,3),rc=rc)
         call LVT_verify(rc,'please specify '//trim(segment_name)//' last grid point lat:')
      enddo
      
      call ESMF_ConfigFindLabel(LVT_config,trim(segment_name)//" last grid point lon:",rc=rc)
      do i=1,LVT_rc%nnest
         call ESMF_ConfigGetAttribute(LVT_config,domain_info(i,4),rc=rc)
         call LVT_verify(rc,'please specify '//trim(segment_name)//' last grid point lon:')
      enddo
      
      call ESMF_ConfigFindLabel(LVT_config,trim(segment_name)//" resolution dlon:",rc=rc)
      do i=1,LVT_rc%nnest
         call ESMF_ConfigGetAttribute(LVT_config,domain_info(i,5),rc=rc)
         call LVT_verify(rc,'please specify '//trim(segment_name)//' resolution dlon:')
      enddo
      
      call ESMF_ConfigFindLabel(LVT_config,trim(segment_name)//" number of lat circles:",rc=rc)
      do i=1,LVT_rc%nnest
         call ESMF_ConfigGetAttribute(LVT_config,domain_info(i,6),rc=rc)
         call LVT_verify(rc,'please specify '//trim(segment_name)//' number of lat circles:')
      enddo

! ---
    case ( "polar" )  
      print*, 'not supported currently'
      stop
#if 0       
      call ESMF_ConfigFindLabel(LVT_config,trim(segment_name)//" lower left lat:",rc=rc)
      do i=1,LVT_rc%nnest
         call ESMF_ConfigGetAttribute(LVT_config,domain_info(i,1),rc=rc)
         call LVT_verify(rc,'please specify '//trim(segment_name)//' lower left lat:')
      enddo
      
      call ESMF_ConfigFindLabel(LVT_config,trim(segment_name)//" lower left lon:",rc=rc)
      do i=1,LVT_rc%nnest
         call ESMF_ConfigGetAttribute(LVT_config,domain_info(i,2),rc=rc)
         call LVT_verify(rc,'please specify '//trim(segment_name)//' lower left lon:')
      enddo
      call ESMF_ConfigFindLabel(LVT_config,trim(segment_name)//" true lat:",rc=rc)
      do i=1,LVT_rc%nnest
         call ESMF_ConfigGetAttribute(LVT_config,domain_info(i,3),rc=rc)
         call LVT_verify(rc,'please specify '//trim(segment_name)//' true lat:')
      enddo
      
      call ESMF_ConfigFindLabel(LVT_config,trim(segment_name)//" standard lon:",rc=rc)
      do i=1,LVT_rc%nnest
         call ESMF_ConfigGetAttribute(LVT_config,domain_info(i,4),rc=rc)
         call LVT_verify(rc,'please specify '//trim(segment_name)//' standard lon:')
      enddo
      
      call ESMF_ConfigFindLabel(LVT_config,trim(segment_name)//" orientation:",rc=rc)
      do i=1,LVT_rc%nnest
         call ESMF_ConfigGetAttribute(LVT_config,domain_info(i,5),rc=rc)
         call LVT_verify(rc,'please specify '//trim(segment_name)//' orientation:')
      enddo
      
      call ESMF_ConfigFindLabel(LVT_config,trim(segment_name)//" resolution:",rc=rc)
      do i=1,LVT_rc%nnest
         call ESMF_ConfigGetAttribute(LVT_config,domain_info(i,6),rc=rc)
         call LVT_verify(rc,'please specify '//trim(segment_name)//' resolution:')
      enddo
      
      call ESMF_ConfigFindLabel(LVT_config,trim(segment_name)//" x-dimension size:",rc=rc)
      do i=1,LVT_rc%nnest
         call ESMF_ConfigGetAttribute(LVT_config,domain_info(i,7),rc=rc)
         call LVT_verify(rc,'please specify '//trim(segment_name)//' x-dimension size:')
      enddo

      call ESMF_ConfigFindLabel(LVT_config,trim(segment_name)//" y-dimension size:",rc=rc)
      do i=1,LVT_rc%nnest
         call ESMF_ConfigGetAttribute(LVT_config,domain_info(i,8),rc=rc)
         call LVT_verify(rc,'please specify '//trim(segment_name)//' y-dimension size:')
      enddo
#endif
! ---
    case ( "UTM" )
      
      call ESMF_ConfigFindLabel(LVT_config,trim(segment_name)//" UTM zone:",rc=rc)
      do i=1,LVT_rc%nnest
         call ESMF_ConfigGetAttribute(LVT_config,domain_info(i,1),rc=rc)
         call LVT_verify(rc,trim(segment_name)//' UTM zone: not defined')
      enddo
      call ESMF_ConfigFindLabel(LVT_config,trim(segment_name)//" northing of SW corner:",rc=rc)
      do i=1,LVT_rc%nnest
         call ESMF_ConfigGetAttribute(LVT_config,domain_info(i,2),rc=rc)
         call LVT_verify(rc,trim(segment_name)//' northing of SW corner: not defined')
      enddo
      
      call ESMF_ConfigFindLabel(LVT_config,trim(segment_name)//" easting of SW corner:",rc=rc)
      do i=1,LVT_rc%nnest
         call ESMF_ConfigGetAttribute(LVT_config,domain_info(i,3),rc=rc)
         call LVT_verify(rc,trim(segment_name)//' easting of SW corner: not defined')
      enddo
      
      call ESMF_ConfigFindLabel(LVT_config,trim(segment_name)//" x-dimension size:",rc=rc)
      do i=1,LVT_rc%nnest
         call ESMF_ConfigGetAttribute(LVT_config,domain_info(i,4),rc=rc)
         call LVT_verify(rc,trim(segment_name)//' x-dimension size: not defined')
      enddo
      
      call ESMF_ConfigFindLabel(LVT_config,trim(segment_name)//" y-dimension size:",rc=rc)
      do i=1,LVT_rc%nnest
         call ESMF_ConfigGetAttribute(LVT_config,domain_info(i,5),rc=rc)
         call LVT_verify(rc,trim(segment_name)//' y-dimension size: not defined')
      enddo
      
      call ESMF_ConfigFindLabel(LVT_config,trim(segment_name)//" resolution:",rc=rc)
      do i=1,LVT_rc%nnest
         call ESMF_ConfigGetAttribute(LVT_config,domain_info(i,6),rc=rc)
         call LVT_verify(rc,trim(segment_name)//' resolution: not defined')
      enddo

!- Non-supported projections (at this time):
    case default 
      print*, 'Reading configuration settings for this projection, ',param_proj,','
      print*, 'is not supported'
      print*, 'Program stopping ....'
      stop
  end select

 end subroutine readDomainConfigSpecs

!BOP
!
! !ROUTINE: read2Ddata
!  \label{read2Ddata}
!
! !INTERFACE:
 subroutine read2Ddata(n, ftn, param_proj, gridDesc, array)

! !USES:
   use map_utils

   implicit none
! !ARGUMENTS:
   integer,  intent(IN)      :: n
   integer,  intent(IN)      :: ftn
   character(50),intent(IN)  :: param_proj
   real,     intent(IN)      :: gridDesc(20)
   real,     intent(INOUT)   :: array(LVT_rc%lnc(n),LVT_rc%lnr(n))
!
! !DESCRIPTION:
!   This routine retrieves a 2-d data from a binary direct access file.
!EOP
   real           :: rlat(LVT_rc%lnc(n),LVT_rc%lnr(n))
   real           :: rlon(LVT_rc%lnc(n),LVT_rc%lnr(n))
   integer*8      :: c, r
   integer*8      :: nc_dom
   integer*8      :: glnc, glnr
   integer*8      :: line1, line2, line
   integer        :: istat
   character*100  :: message  (20)
   real           :: ctmp, rtmp
   type(proj_info) :: proj
! _____________________________________________________________________

! ---  Parameter Projection  ---
  select case ( param_proj )

!- Lat/lon:
   case ( "latlon" )

      do r=1,LVT_rc%lnr(n)
         do c=1,LVT_rc%lnc(n)
            call ij_to_latlon(LVT_domain(n)%lvtproj,float(c),float(r),&
                 rlat(c,r),rlon(c,r))
         enddo
      enddo

      nc_dom = nint((gridDesc(8)-gridDesc(5))/(gridDesc(9)))+1
      do r=1,LVT_rc%lnr(n)
         do c=1,LVT_rc%lnc(n)
            line1 = nint((rlat(c,r)-gridDesc(4))/gridDesc(10))+1 !lat line
            line2 = nint((rlon(c,r)-gridDesc(5))/gridDesc(9))+1  !lon line
            line = (line1-1)*nc_dom + line2
            read(ftn,rec=line) array(c,r)
         enddo
      enddo

!      nc_dom = nint((gridDesc(4)-gridDesc(2))/(gridDesc(5)))+1
!      do r=1,LVT_rc%lnr(n)
!         do c=1,LVT_rc%lnc(n)
!            line1 = nint((rlat(c,r)-gridDesc(1))/gridDesc(6))+1
!            line2 = nint((rlon(c,r)-gridDesc(2))/gridDesc(5))+1
!            line = (line1-1)*nc_dom + line2
!            read(ftn,rec=line) array(c,r)
!         enddo
!      enddo

#if 0
!- Gaussian:
   case ( "gaussian" )
      line1 = gaussian_find_row(LVT_rc%gridDesc(n,4))  -   &
           gaussian_find_row(LVT_rc%gridDesc(n,44)) + 1
      line2 = gaussian_find_col(LVT_rc%gridDesc(n,5))  -   &
           gaussian_find_col(LVT_rc%gridDesc(n,45)) + 1

      do r=1,LVT_rc%lnr(n)
         do c=1,LVT_rc%lnc(n)
            glnc = line2+c-1
            glnr = line1+r-1
            line = (glnr-1)*nint(LVT_rc%gridDesc(n,42))+glnc
            read(ftn,rec=line,iostat=istat) array(c,r)
            if( istat .ne. 0 ) then
               message(1) = 'program:  LIS'
               message(2) = '  routine:  read2DData'
               message(3) = '  iostat != 0'
               call LVT_abort( message )
               call LVT_endrun
            endif
         enddo
      enddo

!- Polar steregraphic:
   case ( "polar" )

      call map_set(PROJ_PS,LVT_rc%gridDesc(n,4),LVT_rc%gridDesc(n,5), &
                   LVT_rc%gridDesc(n,8)*1000.0,                       &
                   LVT_rc%gridDesc(n,11),LVT_rc%gridDesc(n,10),0.0,   &
                   LVT_rc%lnc(n),LVT_rc%lnr(n),proj)

      do r=1,LVT_rc%lnr(n)
         do c=1,LVT_rc%lnc(n)
            call ij_to_latlon(proj,float(c),float(r),rlat(c,r),rlon(c,r))
         enddo
      enddo

      call map_set(PROJ_PS,LVT_rc%lc_gridDesc(n,1),LVT_rc%lc_gridDesc(n,2),  &
                   LVT_rc%lc_gridDesc(n,6)*1000.0,                           &
                   LVT_rc%lc_gridDesc(n,4),LVT_rc%lc_gridDesc(n,3),0.0,      &
                   int(LVT_rc%lc_gridDesc(n,7)),int(LVT_rc%lc_gridDesc(n,8)),&
                   proj)

      do r=1,LVT_rc%lnr(n)
         do c=1,LVT_rc%lnc(n)
            call latlon_to_ij(proj,rlat(c,r),rlon(c,r),ctmp,rtmp)

            line1 = nint(rtmp)
            line2 = nint(ctmp)
            line = (line1-1)*LVT_rc%lc_gridDesc(n,7)+line2
            read(ftn,rec=line) array(c,r)
         enddo
      enddo

   case ( "UTM" )
!rlat/rlon used here to store northing and easting
      do r=1,LVT_rc%lnr(n)
         do c=1,LVT_rc%lnc(n)
            rlat(c,r) = LVT_rc%gridDesc(n,4)+(r-1)*LVT_rc%gridDesc(n,9)
            rlon(c,r) = LVT_rc%gridDesc(n,5)+(c-1)*LVT_rc%gridDesc(n,9)
         enddo
      enddo
      nc_dom = gridDesc(4)

      do r=1,LVT_rc%lnr(n)
         do c=1,LVT_rc%lnc(n)
            line1 = nint((rlat(c,r)-gridDesc(2))/gridDesc(6))+1
            line2 = nint((rlon(c,r)-gridDesc(3))/gridDesc(6))+1
            line = (line1-1)*nc_dom +line2
           read(ftn,rec=line) array(c,r)
         enddo
      enddo
#endif
   case default
      write(*,*) 'This parameter projection is not supported...'
      write(*,*) 'Program stopping ....'
      stop

   end select

 end subroutine read2Ddata


!BOP
! !ROUTINE: readLCdata
!  \label{rreadLCdata}
! 
! !INTERFACE: 
 subroutine readLCdata( n, ftn, param_proj, param_gridDesc, array, domain )

! !USES: 

   implicit none

! !ARGUMENTS: 
   integer,  intent(IN)     :: n                  ! Nest
   integer,  intent(IN)     :: ftn                ! File number
   character(50),intent(IN) :: param_proj         ! LC file projection
   real,     intent(IN)     :: param_gridDesc(20)  ! LC file grid description array
   real,     intent(INOUT)  :: array(LVT_rc%lnc(n),LVT_rc%lnr(n),LVT_rc%nt)  ! Read in LC array
   integer,  intent(IN), optional :: domain       ! Mask domain information - ?
!
! !DESCRIPTION: 
!  This routine reads a 3-d vegetation distribution data from a binary 
!   sequential access file. 
!EOP  
   real         :: rlat(LVT_rc%lnc(n),LVT_rc%lnr(n))
   real         :: rlon(LVT_rc%lnc(n),LVT_rc%lnr(n))
   integer      :: c, r, t
   integer      :: line1, line2, line
   integer      :: glnc, glnr, num_lon_pts
   integer      :: pnc, pnr
   real, allocatable :: tsum(:,:)
   real, allocatable :: vegtype(:,:)
! _________________________________________________________________________________

   pnr = nint((param_gridDesc(7)-param_gridDesc(4))/param_gridDesc(10)) + 1
   pnc = nint((param_gridDesc(8)-param_gridDesc(5))/param_gridDesc(9)) + 1

!- Lat/lon projection:
   if( param_proj == "latlon" ) then 

      do r=1,LVT_rc%lnr(n)
         do c=1,LVT_rc%lnc(n)
            call ij_to_latlon(LVT_domain(n)%lvtproj,float(c),float(r),&
                 rlat(c,r),rlon(c,r))
         enddo
      enddo
      
   !- 3D - tiled:
      if( param_gridDesc(9).ne.0.01 ) then
#if defined ( INC_WATER_PTS )
       do t=1,LVT_rc%nt-1
#else
       do t=1,LVT_rc%nt
#endif
         do r=1,LVT_rc%lnr(n)
            do c=1,LVT_rc%lnc(n)
               line1 = nint((rlat(c,r)-param_gridDesc(4))/param_gridDesc(10))+1
               line2 = nint((rlon(c,r)-param_gridDesc(5))/param_gridDesc(9))+1
               line = (line1-1)*pnc + line2 + (t-1)*(pnc*pnr)

            !- Read in 3D parameter data:
               read(ftn,rec=line) array(c,r,t)
            enddo
         enddo
       enddo

   ! - 2D - 1KM:
      else
         allocate(tsum(LVT_rc%lnc(n),LVT_rc%lnr(n)))
         do r=1,LVT_rc%lnr(n)
            do c=1,LVT_rc%lnc(n)
               line1 = nint((rlat(c,r)-param_gridDesc(4))/param_gridDesc(10))+1
               line2 = nint((rlon(c,r)-param_gridDesc(5))/param_gridDesc(9))+1
               line = (line1-1)*pnc + line2
               read(ftn,rec=line) tsum(c,r)
#if ( defined INC_WATER_PTS )
               !kludged to include water points
               ! The right way to do this would be read an appropriate mask and veg
               ! file with water points
               if(tsum(c,r).le.0) then
                  tsum(c,r) = LVT_rc%waterclass
               endif
               !kludge
#endif
               if(nint(tsum(c,r)).ne.LVT_rc%waterclass) then
                  array(c,r,NINT(tsum(c,r))) = 1.0
               endif
            enddo
         enddo
         deallocate(tsum)

      end if

! *********************************************************************
#if 0
!- Gaussian:
   elseif( param_proj == "gaussian" ) then 

      glnr = gaussian_find_row(LVT_rc%gridDesc(n,4)) -   &
              gaussian_find_row(LVT_rc%gridDesc(n,44)) + 1
      glnc = gaussian_find_col(LVT_rc%gridDesc(n,5)) -   &
              gaussian_find_col(LVT_rc%gridDesc(n,45)) + 1
      num_lon_pts = nint(360/(param_gridDesc(5)))

      do t=1,LVT_rc%nt
         do r=1,LVT_rc%lnr(n)
            do c=1,LVT_rc%lnc(n)
               line1 = glnr+r-1
               line2 = glnc+c-1
               line = (line1-1)*num_lon_pts+line2+(t-1)*(pnc*pnr)

            !- Read in 3D parameter data:
               read(ftn,rec=line) array(c,r,t)
            enddo
         enddo
      enddo

! *********************************************************************

!- UTM:
   elseif( param_proj == "UTM" ) then

      allocate(vegtype(LVT_rc%lnc(n), LVT_rc%lnr(n)))

      pnc = param_gridDesc(4)

    !-rlat/rlon used here to store northing and easting
      do r=1,LVT_rc%lnr(n)
         do c=1,LVT_rc%lnc(n)
            rlat(c,r) = LVT_rc%gridDesc(n,4)+(r-1)*LVT_rc%gridDesc(n,9)
            rlon(c,r) = LVT_rc%gridDesc(n,5)+(c-1)*LVT_rc%gridDesc(n,9)
         enddo
      enddo

      do r=1,LVT_rc%lnr(n)
         do c=1,LVT_rc%lnc(n)
            line1 = nint((rlat(c,r)-param_gridDesc(2))/param_gridDesc(6))+1
            line2 = nint((rlon(c,r)-param_gridDesc(3))/param_gridDesc(6))+1
            line = (line1-1)*pnc +line2

         !- Read in 2D parameter data:
            read(ftn,rec=line) vegtype(c,r)

            if(nint(vegtype(c,r)).ne.LVT_rc%waterclass.and.&
               nint(vegtype(c,r)).ne.LVT_rc%udef) then
               array(c,r,NINT(vegtype(c,r))) = 1.0
            endif
         enddo
      enddo
      deallocate(vegtype)

#endif

   else 
      write(*,*) 'This parameter projection is not supported...'
      write(*,*) 'Program stopping ....'
      stop
   endif

 end subroutine readLCdata

!BOP
! !ROUTINE: paramProcWrite
! \label{paramProcWrite}
!
! !INTERFACE: 
  subroutine paramProcWrite()

! !USES:
    integer    :: n 

    integer               :: ftn
    integer               :: iret
    integer               :: m
    integer               :: dimID(3),monthID,qID
    integer               :: tdimID, xtimeID
    character(len=8)      :: date
    character(len=10)     :: time
    character(len=5)      :: zone
    integer, dimension(8) :: values
    integer               :: varid
    integer               :: c,r
    integer               :: shuffle, deflate, deflate_level
    integer               :: xlatid, xlonid
    type(LVT_paramEntry)  :: xlat, xlon
! ________________________________________________________

    n = 1
    shuffle = 1
    deflate = 1
    deflate_level =9

    iret=nf90_create(path=trim(LVT_LSMparam_struc(n)%param_filename),&
         cmode=nf90_clobber, ncid=ftn)
    call LVT_verify(iret,'creating netcdf file failed in LVT_LSMparamProcWrite')

!-- General header:
    call date_and_time(date,time,zone,values)

!-- Write out dimensions headers:

 !- Grid-domain dimensions:
    call verify(nf90_def_dim(ftn,'east_west',LVT_rc%gnc(n),dimID(1)))
    call verify(nf90_def_dim(ftn,'north_south',LVT_rc%gnr(n),dimID(2)))
    
 !- Time-based dimensions:
    call verify(nf90_def_dim(ftn,'time',1,tdimID))
    call verify(nf90_def_var(ftn,'time',nf90_float,dimids=tdimID,&
         varID=xtimeID))

 !- LIS-Domain Grid file attributes:
    select case ( LVT_rc%lis_map_proj )

     case ( "latlon" )
       call verify(nf90_put_att(ftn,NF90_GLOBAL,"MAP_PROJECTION", &
            "EQUIDISTANT CYLINDRICAL"))
       call verify(nf90_put_att(ftn,NF90_GLOBAL,"SOUTH_WEST_CORNER_LAT", &
            LVT_rc%gridDesc(n,4)))
       call verify(nf90_put_att(ftn,NF90_GLOBAL,"SOUTH_WEST_CORNER_LON", &
            LVT_rc%gridDesc(n,5)))
       call verify(nf90_put_att(ftn,NF90_GLOBAL,"DX", &
            LVT_rc%gridDesc(n,9)))
       call verify(nf90_put_att(ftn,NF90_GLOBAL,"DY", &
            LVT_rc%gridDesc(n,10)))       
     case ( "mercator" )
       call verify(nf90_put_att(ftn,NF90_GLOBAL,"MAP_PROJECTION", &
            "MERCATOR"))
       call verify(nf90_put_att(ftn,NF90_GLOBAL,"SOUTH_WEST_CORNER_LAT", &
            LVT_rc%gridDesc(n,4)))
       call verify(nf90_put_att(ftn,NF90_GLOBAL,"SOUTH_WEST_CORNER_LON", &
            LVT_rc%gridDesc(n,5)))
       call verify(nf90_put_att(ftn,NF90_GLOBAL,"TRUELAT1", &
            LVT_rc%gridDesc(n,10)))
       call verify(nf90_put_att(ftn,NF90_GLOBAL,"STANDARD_LON", &
            LVT_rc%gridDesc(n,11)))
       call verify(nf90_put_att(ftn,NF90_GLOBAL,"DX", &
            LVT_rc%gridDesc(n,8)))
       call verify(nf90_put_att(ftn,NF90_GLOBAL,"DY", &
            LVT_rc%gridDesc(n,9)))
     case ( "lambert" )    ! Lambert conformal
       call verify(nf90_put_att(ftn,NF90_GLOBAL,"MAP_PROJECTION", &
            "LAMBERT CONFORMAL"))
       call verify(nf90_put_att(ftn,NF90_GLOBAL,"SOUTH_WEST_CORNER_LAT", &
            LVT_rc%gridDesc(n,4)))
       call verify(nf90_put_att(ftn,NF90_GLOBAL,"SOUTH_WEST_CORNER_LON", &
            LVT_rc%gridDesc(n,5)))
       call verify(nf90_put_att(ftn,NF90_GLOBAL,"TRUELAT1", &
            LVT_rc%gridDesc(n,10)))
       call verify(nf90_put_att(ftn,NF90_GLOBAL,"TRUELAT2", &
            LVT_rc%gridDesc(n,7)))
       call verify(nf90_put_att(ftn,NF90_GLOBAL,"STANDARD_LON", &
            LVT_rc%gridDesc(n,11)))
       call verify(nf90_put_att(ftn,NF90_GLOBAL,"DX", &
            LVT_rc%gridDesc(n,8)))
       call verify(nf90_put_att(ftn,NF90_GLOBAL,"DY", &
            LVT_rc%gridDesc(n,9)))
     case ( "polar" )    ! Polar stereographic
       call verify(nf90_put_att(ftn,NF90_GLOBAL,"MAP_PROJECTION", &
            "POLAR STEREOGRAPHIC"))
       call verify(nf90_put_att(ftn,NF90_GLOBAL,"SOUTH_WEST_CORNER_LAT", &
            LVT_rc%gridDesc(n,4)))
       call verify(nf90_put_att(ftn,NF90_GLOBAL,"SOUTH_WEST_CORNER_LON", &
            LVT_rc%gridDesc(n,5)))
       call verify(nf90_put_att(ftn,NF90_GLOBAL,"TRUELAT1", &
            LVT_rc%gridDesc(n,10)))
       call verify(nf90_put_att(ftn,NF90_GLOBAL,"ORIENT", &
            LVT_rc%gridDesc(n,7)))
       call verify(nf90_put_att(ftn,NF90_GLOBAL,"STANDARD_LON", &
            LVT_rc%gridDesc(n,11)))
       call verify(nf90_put_att(ftn,NF90_GLOBAL,"DX", &
            LVT_rc%gridDesc(n,8)))
       call verify(nf90_put_att(ftn,NF90_GLOBAL,"DY", &
            LVT_rc%gridDesc(n,9)))
    end select

 !- Include water points:
#if ( defined INC_WATER_PTS )
    call verify(nf90_put_att(ftn,NF90_GLOBAL,"INC_WATER_PTS", &
        "true"))
#else
    call verify(nf90_put_att(ftn,NF90_GLOBAL,"INC_WATER_PTS", &
        "false"))
# endif

! - Write Parameter Header Data:
    call writeParamHeaders(n, ftn, dimID, monthID, qID)

! - write lat, lons
    xlat%short_name    = "lat"
    xlat%standard_name = "latitude"
    xlat%units         = "degrees_north"
    allocate(xlat%value(LVT_rc%lnc(n),LVT_rc%lnr(n),1))

    xlon%short_name    = "lon"
    xlon%standard_name = "longitude"
    xlon%units         = "degrees_east"
    allocate(xlon%value(LVT_rc%lnc(n),LVT_rc%lnr(n),1))

    do r=1,LVT_rc%lnr(n)
       do c=1,LVT_rc%lnc(n)
          call ij_to_latlon(LVT_domain(n)%lvtproj,&
               real(c), real(r), xlat%value(c,r,1),&
               xlon%value(c,r,1))
       enddo
    enddo
    
    call LVT_verify(nf90_def_var(ftn,trim(xlat%short_name),&
         nf90_float, dimids=dimID(1:2),varid=xlatid),&
         'nf90_def_var failed for xlat')

    call LVT_verify(nf90_put_att(ftn,xlatid, &
         "standard_name",trim(xlat%standard_name)),&
         'nf90_put_att failed for xlat:standard_name')
    call LVT_verify(nf90_put_att(ftn,xlatid, &
         "units",trim(xlat%units)),&
         'nf90_put_att failed for xlat:units')
    call LVT_verify(nf90_put_att(ftn,xlatid, &
         "scale_factor",1.0),&
         'nf90_put_att failed for xlat:scale_factor')
    call LVT_verify(nf90_put_att(ftn,xlatid, &
         "add_offset",0.0),&
         'nf90_put_att failed for xlat:add_offset')
    call LVT_verify(nf90_put_att(ftn,xlatid, &
         "missing_value",LVT_rc%udef),&
         'nf90_put_att failed for xlat:missing_value')
    call LVT_verify(nf90_put_att(ftn,xlatid, &
         "vmin",0.0),&
         'nf90_put_att failed for xlat:vmin')
    call LVT_verify(nf90_put_att(ftn,xlatid, &
         "vmax",0.0),&
         'nf90_put_att failed for xlat:vmax')

    call LVT_verify(nf90_def_var(ftn,trim(xlon%short_name),&
         nf90_float, dimids=dimID(1:2),varid=xlonid),&
         'nf90_def_var failed for xlon')

    call LVT_verify(nf90_put_att(ftn,xlonid, &
         "standard_name",trim(xlon%standard_name)),&
         'nf90_put_att failed for xlon:standard_name')
    call LVT_verify(nf90_put_att(ftn,xlonid, &
         "units",trim(xlon%units)),&
         'nf90_put_att failed for xlon:units')
    call LVT_verify(nf90_put_att(ftn,xlonid, &
         "scale_factor",1.0),&
         'nf90_put_att failed for xlon:scale_factor')
    call LVT_verify(nf90_put_att(ftn,xlonid, &
         "add_offset",0.0),&
         'nf90_put_att failed for xlon:add_offset')
    call LVT_verify(nf90_put_att(ftn,xlonid, &
         "missing_value",LVT_rc%udef),&
         'nf90_put_att failed for xlon:missing_value')
    call LVT_verify(nf90_put_att(ftn,xlonid, &
         "vmin",0.0),&
         'nf90_put_att failed for xlon:vmin')
    call LVT_verify(nf90_put_att(ftn,xlonid, &
         "vmax",0.0),&
         'nf90_put_att failed for xlon:vmax')

    ! - Write file attributes information:
    call verify(nf90_put_att(ftn,NF90_GLOBAL,"title", &
         "Land Verification Toolkit (LVT) output"))
    call verify(nf90_put_att(ftn,NF90_GLOBAL,"institution", &
         "NASA GSFC Hydrological Sciences Laboratory"))
    call verify(nf90_put_att(ftn,NF90_GLOBAL,"history", &
         "created on date: "//date(1:4)//"-"//date(5:6)//"-"//&
         date(7:8)//"T"//time(1:2)//":"//time(3:4)//":"//time(5:10)))
    call verify(nf90_put_att(ftn,NF90_GLOBAL,"references", &
         "Kumar_etal_EMS_2006, Kumar_etal_GMD_2012"))
    call verify(nf90_put_att(ftn,NF90_GLOBAL,"comment", &
         "website: http://lis.gsfc.nasa.gov/LVT"))

    call verify(nf90_enddef(ftn))

    call LVT_verify(nf90_put_var(ftn,xlatid, xlat%value(:,:,1),&
         (/1,1/),(/LVT_rc%gnc(n),LVT_rc%gnr(n)/)),&
         'nf90_put_att failed for xlat')

    call LVT_verify(nf90_put_var(ftn,xlonid, xlon%value(:,:,1),&
         (/1,1/),(/LVT_rc%gnc(n),LVT_rc%gnr(n)/)),&
         'nf90_put_att failed for xlon')

    call verify(nf90_put_var(ftn,xtimeID,0.0))

! - Write Parameter Output Data:
    call writeParamData(n,ftn)

! - Close file:
    call verify(nf90_close(ftn))


    write(*,*) "-- Wrote parameters to netcdf output file -- "

  end subroutine paramProcWrite

  subroutine writeParamHeaders(n, ftn, dimID, monthID, qID)
    
! !USES:
    integer     :: n 
    integer     :: ftn
    integer     :: dimID(3)
    integer     :: monthID
    integer     :: qID
    integer     :: k
    
    call LMLC_writeHeader(n,ftn,dimID)

  end subroutine writeParamHeaders

  subroutine writeParamData(n, ftn)

! !USES:

    integer  :: n 
    integer  :: ftn
    integer  :: k

    call LMLC_writeData(n,ftn)

  end subroutine writeParamData

  subroutine LMLC_writeHeader(n,ftn,dimID)
    

    integer      :: n 
    integer      :: ftn
    integer      :: dimID(3)
    integer      :: tdimID(3)

    tdimID(1) = dimID(1)
    tdimID(2) = dimID(2)
   
    call verify(nf90_def_dim(ftn,'sfctypes',&
         LVT_LSMparam_struc(n)%landcover%vlevels,tdimID(3)))

    call LVT_writeNETCDFdataHeader(n,ftn,tdimID,&
         LVT_LSMparam_struc(n)%landmask)
    call LVT_writeNETCDFdataHeader(n,ftn,tdimID,&
         LVT_LSMparam_struc(n)%landcover)
    call LVT_writeNETCDFdataHeader(n,ftn,tdimID,&
         LVT_LSMparam_struc(n)%sfctype)

    call verify(nf90_put_att(ftn,NF90_GLOBAL,"BARESOILCLASS", &
         LVT_rc%bareclass))
    call verify(nf90_put_att(ftn,NF90_GLOBAL,"URBANCLASS", &
         LVT_rc%urbanclass))
    call verify(nf90_put_att(ftn,NF90_GLOBAL,"SNOWCLASS", &
         LVT_rc%snowclass))
    call verify(nf90_put_att(ftn,NF90_GLOBAL,"WATERCLASS", &
         LVT_rc%waterclass))
    call verify(nf90_put_att(ftn,NF90_GLOBAL,"WETLANDCLASS", &
         LVT_rc%wetlandclass))
    call verify(nf90_put_att(ftn,NF90_GLOBAL,"GLACIERCLASS", &
         LVT_rc%glacierclass))

! - Enter number of vegetation land use types only:
    call verify(nf90_put_att(ftn,NF90_GLOBAL,"NUMVEGTYPES", &
         LVT_rc%nt))

  end subroutine LMLC_writeHeader

  subroutine LMLC_writeData(n,ftn)


    integer  :: n 
    integer  :: ftn

    call LVT_writeNETCDFdata(n,ftn,LVT_LSMparam_struc(n)%landmask)
    call LVT_writeNETCDFdata(n,ftn,LVT_LSMparam_struc(n)%landcover)
    call LVT_writeNETCDFdata(n,ftn,LVT_LSMparam_struc(n)%sfctype)

  end subroutine LMLC_writeData

  subroutine LVT_writeNETCDFdataHeader(n,ftn,dimID, paramEntry)
    
    integer               :: n 
    integer               :: ftn
    integer               :: dimID(3)
    integer               :: shuffle, deflate, deflate_level
    type(LVT_paramEntry)        :: paramEntry
    
    shuffle = 1
    deflate = 1
    deflate_level =9

    if(paramEntry%selectOpt.gt.0) then 
       if(paramEntry%vlevels.gt.1) then 
          call verify(nf90_def_var(ftn,trim(paramEntry%short_name),&
               nf90_float, dimids = dimID, varID=paramEntry%vid))
       else
          call verify(nf90_def_var(ftn,trim(paramEntry%short_name),&
               nf90_float, dimids = dimID(1:2), varID=paramEntry%vid))
       endif

    endif
  end subroutine LVT_writeNETCDFdataHeader

  subroutine LVT_writeNETCDFdata( n, ftn, paramEntry )
    
   integer              :: n 
   integer              :: ftn
   type(LVT_paramEntry) :: paramEntry    

   if(paramEntry%selectOpt.gt.0) then 
      
      if(paramEntry%vlevels.gt.1) then 
         call verify(nf90_put_var(ftn, paramEntry%vid,paramEntry%value,&
              (/1,1,1/),(/LVT_rc%gnc(n),LVT_rc%gnr(n),&
              paramEntry%vlevels/)))
      else
         call LVT_verify(nf90_put_var(ftn, paramEntry%vid,paramEntry%value,&
              (/1,1/),(/LVT_rc%gnc(n),LVT_rc%gnr(n)/)),&
              'nf90_put_var failed in LVT_writeNETCDFdata')
      endif

#if 0 
      call verify(nf90_put_att(ftn,paramEntry%vid,&
           "standard_name",trim(paramEntry%standard_name)))
      call verify(nf90_put_att(ftn,paramEntry%vid,&
           "units",trim(paramEntry%units)))
      call verify(nf90_put_att(ftn,paramEntry%vid,&
           "scale_factor",1.0))
      call verify(nf90_put_att(ftn,paramEntry%vid,&
           "add_offset",0.0))
      call verify(nf90_put_att(ftn,paramEntry%vid,&
           "missing_value",LVT_rc%udef))
!      call verify(nf90_put_att(ftn,paramEntry%vid,&
!           "_FillValue",LVT_rc%udef))
      call verify(nf90_put_att(ftn,paramEntry%vid,&
           "vmin",paramEntry%valid_min))
      call verify(nf90_put_att(ftn,paramEntry%vid,&
           "vmax",paramEntry%valid_max))
      call verify(nf90_put_att(ftn,paramEntry%vid,&
           "num_bins",paramEntry%num_bins))
#endif
   endif

 end subroutine LVT_writeNETCDFdata

!BOP
! !ROUTINE: verify
! 
! !INTERFACE:
  subroutine verify(ierr)
! 
! !ARGUMENTS:
    implicit none
    integer, intent(in)          :: ierr
! 
! !DESCRIPTION:
! This is an error check routine. Program exits in case of error
! with an associated error message written to the 'abort message'
! file. 
! 
!EOP

    if ( ierr /= 0 ) then
       print*,'ERR: in LVT: Stopping.'
       stop
    endif

  end subroutine Verify


!BOP
! !ROUTINE: verify
! 
! !INTERFACE:
  subroutine LVT_verify(ierr,msg)
! 
! !ARGUMENTS:
    implicit none
    integer, intent(in)          :: ierr
    character(len=*), intent(in) :: msg  
! 
! !DESCRIPTION:
! This is an error check routine. Program exits in case of error
! with an associated error message written to the 'abort message'
! file. 
! 
!EOP

    if ( ierr /= 0 ) then
       print*,'ERR: ',msg,' Stopping.'
       stop
    endif

  end subroutine LVT_verify
  
end module preprocMod
