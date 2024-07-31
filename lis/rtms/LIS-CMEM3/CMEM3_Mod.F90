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
module CMEM3_Mod
!BOP
!
! !MODULE: CMEM3_Mod
!
! !DESCRIPTION:
!    This module provides the routines to control the execution of 
!    LIS-CMEM3 from within LIS. The subroutines provide mapping of the surface
!    and atmospheric profiles to CMEM3, which then runs the forward 
!    model to generate the radiative quantities. 
!
! !REVISION HISTORY:
! 03 Feb 2011: Yudong Tian; Modifed from CRTM2_handlerMod.F90 to support CMEM3 
! 23 Feb 2011: Yudong Tian; Added veg_frac in cmem_struc to save it in output 
!
! !USES:        


#if (defined RTMS)

  use ESMF
  use LIS_constantsMod, only : LIS_CONST_PATH_LEN

  implicit none
 
  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: CMEM3_initialize
  public :: CMEM3_f2t
  public :: CMEM3_run
  public :: CMEM3_geometry 
!!$  public :: CMEM3_landmatch
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  public :: cmem3_struc
!EOP
  type, public ::  cmem3_type_dec 
     character*256                         :: sensor_id
     character(len=LIS_CONST_PATH_LEN)     :: freqfile 
     integer                               :: nfs 
     real                                  :: fghz(100)    ! max 100 freqs supported
     real                                  :: theta(100) 
     character*50                          :: chName(100)
     !-------------- atmosphere -------------------------!
     real, allocatable                         :: tsurf(:)
     real, allocatable                         :: tau_atm0(:)  !atm zenith optical depth
     real, allocatable                         :: tb_ad(:)
     real, allocatable                         :: tb_au(:)
     !-------------- soil -------------------------------!
     real, allocatable                         :: tsoil(:)
     real, allocatable                         :: wc(:)
     real, allocatable                         :: sr(:)
     real, allocatable                         :: srmax(:)  ! testing allowing sr to be linked to soil moisture
     real, allocatable                         :: srmax2srmin(:)
     real, allocatable                         :: sand(:)
     real, allocatable                         :: clay(:)
     real, allocatable                         :: wc_fac(:)  ! soil water content
     real, allocatable                         :: relsmc(:)

     !-------------- veg --------------------------------!
     integer, allocatable                      :: vtype(:)
     real, allocatable                         :: tveg(:)
     real, allocatable                         :: h_veg(:)
     real, allocatable                         :: lai(:) 
     real, allocatable                         :: wc_veg(:)
     real, allocatable                         :: rho_veg(:)
     real, allocatable                         :: m_d(:)
     real, allocatable                         :: d_leaf(:) 
!     real, allocatable                         :: veg_frac(:)    ! not used in CMEM for now
     real, allocatable                         :: vwc2lai(:)     ! VWC / LAI ratio 
     real, allocatable                         :: bgf_fixed(:)     ! bare ground fraction
     real, allocatable                         :: bgf_total(:)         ! bare ground fraction total for grid cell
     real, allocatable                         :: bgf_v(:)         ! bare ground fraction total for veg portion of grid cell
     real, allocatable                         :: k_lai2vgf(:)  ! gvf of vegetated portion = 1*exp(-k_lai2vgf*LAI) 
     !-------------- snow -------------------------------!
     real, allocatable                         :: sn_t(:) 
     real, allocatable                         :: sn_moist(:) 
     real, allocatable                         :: sn_density(:)
     real, allocatable                         :: sn_depth(:)
     real, allocatable                         :: sn_gsize(:)
     real, allocatable                         :: sn_depth_bias_fac(:) ! for perturbation
     !-------------- output -----------------------------!
     real, allocatable                         :: em(:, :)   ! ntiles, nfs 
     real, allocatable                         :: emH(:, :)   ! ntiles, nfs 
     real, allocatable                         :: emV(:, :)
     real, allocatable                         :: toaH(:, :) 
     real, allocatable                         :: toaV(:, :) 
     integer          :: fileopen
  end type cmem3_type_dec

  type, public ::  sm_correction_dec
     integer                               :: c_type
     character(len=LIS_CONST_PATH_LEN)     :: src_mean_file
     character(len=LIS_CONST_PATH_LEN)     :: src_sigma_file
     character(len=LIS_CONST_PATH_LEN)     :: dst_mean_file
     character(len=LIS_CONST_PATH_LEN)     :: dst_sigma_file
     real                                  :: gridDesc(8)
     real, allocatable                         :: src_mean(:, :), src_sigma(:, :)
     real, allocatable                         :: dst_mean(:, :), dst_sigma(:, :)
  end type sm_correction_dec

  type(cmem3_type_dec), allocatable :: cmem3_struc(:) 
  type(sm_correction_dec), allocatable :: sm_correction_struc(:)
  SAVE

contains
!BOP
! 
! !ROUTINE: CMEM3_initialize
! \label{CMEM3_initialize}
! 
! !INTERFACE:
  subroutine CMEM3_initialize()
! !USES:
    use LIS_coreMod,    only : LIS_rc, LIS_config
    use LIS_logMod,     only : LIS_logunit, LIS_verify, LIS_getNextUnitNumber, &
                               LIS_releaseUnitNumber
    use LIS_RTMMod,     only : LIS_sfcState,LIS_forwardState
    use LIS_fileIOMod,  only : LIS_readData, LIS_readDomainConfigSpecs

! !DESCRIPTION:        
!
!  This routine creates the datatypes and allocates memory for noah-specific
!  variables. It also invokes the routine to read the runtime specific options
!  for noah from the configuration file. 
! 
!  The routines invoked are: 
!  \begin{description}
!   \item[readCMEM3crd](\ref{readCMEM3crd}) \newline
!    reads the runtime options for CMEM3 EMonly
!  \end{description}
!EOP
    implicit none

    integer :: n 
    integer :: m 
    integer :: j
    integer :: n_Sensors, n_Channels, n_tiles, n_fs, nc, np, rc, ftn
    integer :: status
    character*100 :: temp, fieldname
    real, allocatable :: tmp_gridDesc(:, :)

    allocate(tmp_gridDesc(LIS_rc%nnest, 8))
    allocate( cmem3_struc(LIS_rc%nnest) )
    allocate(sm_correction_struc(LIS_rc%nnest))

!SVK: commenting this out for now
#if 0 
    call ESMF_ConfigFindLabel(LIS_config,&
         "RTM input soil moisture correction:",rc=rc)
    do n=1,LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,&
            sm_correction_struc(n)%c_type,rc=rc)
    end do
    
    call ESMF_ConfigFindLabel(LIS_config,&
         "RTM input soil moisture correction src mean file:",rc=rc)
    do n=1,LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,&
            sm_correction_struc(n)%src_mean_file,rc=rc)
    end do
    
    call ESMF_ConfigFindLabel(LIS_config,&
         "RTM input soil moisture correction src sigma file:",rc=rc)
    do n=1,LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,&
            sm_correction_struc(n)%src_sigma_file,rc=rc)
    end do
    
    call ESMF_ConfigFindLabel(LIS_config,&
         "RTM input soil moisture correction dst mean file:",rc=rc)
    do n=1,LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,&
            sm_correction_struc(n)%dst_mean_file,rc=rc)
    end do
    
    call ESMF_ConfigFindLabel(LIS_config,&
         "RTM input soil moisture correction dst sigma file:",rc=rc)
    do n=1,LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,&
            sm_correction_struc(n)%dst_sigma_file,rc=rc)
    end do
    
    call LIS_readDomainConfigSpecs("RTM input soil moisture correction", &
         tmp_gridDesc)
    
    do n=1,LIS_rc%nnest
       if (sm_correction_struc(n)%c_type .eq. 1) then
          sm_correction_struc(n)%gridDesc(:) = tmp_gridDesc(n, :)
          allocate(sm_correction_struc(n)%src_mean(&
               LIS_rc%lnc(n),LIS_rc%lnr(n)))
          allocate(sm_correction_struc(n)%src_sigma(&
               LIS_rc%lnc(n),LIS_rc%lnr(n)))
          allocate(sm_correction_struc(n)%dst_mean(&
               LIS_rc%lnc(n),LIS_rc%lnr(n)))
          allocate(sm_correction_struc(n)%dst_sigma(&
               LIS_rc%lnc(n),LIS_rc%lnr(n)))
          
          ftn = LIS_getNextUnitNumber()
          open(ftn, file=sm_correction_struc(n)%src_mean_file, &
               access="direct",  &
               form="unformatted", recl=4)
          call LIS_readData(n, ftn, sm_correction_struc(n)%gridDesc, &
               sm_correction_struc(n)%src_mean)
          call LIS_releaseUnitNumber(ftn)
          
          ftn = LIS_getNextUnitNumber()
          open(ftn, file=sm_correction_struc(n)%src_sigma_file, &
               access="direct",  &
               form="unformatted", recl=4)
          call LIS_readData(n, ftn, sm_correction_struc(n)%gridDesc, &
               sm_correction_struc(n)%src_sigma)
          call LIS_releaseUnitNumber(ftn)
          
          ftn = LIS_getNextUnitNumber()
          open(ftn, file=sm_correction_struc(n)%dst_mean_file, &
               access="direct",  &
               form="unformatted", recl=4)
          call LIS_readData(n, ftn, sm_correction_struc(n)%gridDesc,&
               sm_correction_struc(n)%dst_mean)
          call LIS_releaseUnitNumber(ftn)
          
          ftn = LIS_getNextUnitNumber()
          open(ftn, file=sm_correction_struc(n)%dst_sigma_file, &
               access="direct",  &
               form="unformatted", recl=4)
          call LIS_readData(n, ftn, sm_correction_struc(n)%gridDesc, &
               sm_correction_struc(n)%dst_sigma)
          call LIS_releaseUnitNumber(ftn)
          
       end if
    end do
#endif
    call readCMEM3crd()

    do n=1,LIS_rc%nnest
!---------------------------------------------------------------------
! Allocate ESMF State for mapping surface properties
!---------------------------------------------------------------------

       call add_fields_toState(n,LIS_sfcState(n), "Wind Speed")
       call add_fields_toState(n,LIS_sfcState(n), "Land Coverage")
       call add_fields_toState(n,LIS_sfcState(n), "Land Type")
       call add_fields_toState(n,LIS_sfcState(n), "Land Temperature")
       call add_fields_toState(n,LIS_sfcState(n), "Snow Temperature")
       call add_fields_toState(n,LIS_sfcState(n), "Soil Moisture Content")
       call add_fields_toState(n,LIS_sfcState(n), "Soil Temperature")
       call add_fields_toState(n,LIS_sfcState(n), "Canopy Water Content")
       call add_fields_toState(n,LIS_sfcState(n), "Vegetation Fraction")
       call add_fields_toState(n,LIS_sfcState(n), "Snow Coverage")
       call add_fields_toState(n,LIS_sfcState(n), "Snow Depth")
       call add_fields_toState(n,LIS_sfcState(n), "Snow Density")
       call add_fields_toState(n,LIS_sfcState(n), "Leaf Area Index")
       call add_fields_toState(n,LIS_sfcState(n), "Relative Soil Moisture Content")
       
       do j=1, cmem3_struc(n)%nfs
!for each channel, create a field with both H and V polarizations          
          call add_fields_toState(n,LIS_forwardState(n),&
               trim(cmem3_struc(n)%chName(j))//"_Hpol")
          call add_fields_toState(n,LIS_forwardState(n),&
               trim(cmem3_struc(n)%chName(j))//"_Vpol")
       enddo
       
       ! Number of freqs 
       n_fs = cmem3_struc(n)%nfs
       
       ! Allocate tile space variables 
       n_tiles = LIS_rc%ntiles(n) 
       ! atmos
       allocate( cmem3_struc(n)%tsurf	(n_tiles) ) 
       allocate( cmem3_struc(n)%tau_atm0(n_tiles) ) 
       allocate( cmem3_struc(n)%tb_ad	(n_tiles) ) 
       allocate( cmem3_struc(n)%tb_au	(n_tiles) ) 
       ! soil
       allocate( cmem3_struc(n)%tsoil	(n_tiles) ) 
       allocate( cmem3_struc(n)%wc    	(n_tiles) ) 
       allocate( cmem3_struc(n)%sr   	(n_tiles) ) 
       allocate( cmem3_struc(n)%srmax   (n_tiles) ) 
       allocate( cmem3_struc(n)%srmax2srmin (n_tiles) )  
       allocate( cmem3_struc(n)%sand 	(n_tiles) ) 
       allocate( cmem3_struc(n)%clay 	(n_tiles) ) 
       allocate( cmem3_struc(n)%wc_fac	(n_tiles) ) 
       allocate( cmem3_struc(n)%relsmc	(n_tiles) ) 

       ! veg
       allocate( cmem3_struc(n)%vtype	(n_tiles) ) 
       allocate( cmem3_struc(n)%tveg 	(n_tiles) ) 
       allocate( cmem3_struc(n)%h_veg	(n_tiles) ) 
       allocate( cmem3_struc(n)%lai  	(n_tiles) ) 
       allocate( cmem3_struc(n)%wc_veg  (n_tiles) ) 
       allocate( cmem3_struc(n)%rho_veg (n_tiles) ) 
       allocate( cmem3_struc(n)%m_d  	(n_tiles) ) 
       allocate( cmem3_struc(n)%d_leaf  (n_tiles) ) 
!       allocate( cmem3_struc(n)%veg_frac(n_tiles) ) 
       allocate( cmem3_struc(n)%vwc2lai (n_tiles) ) 
       allocate( cmem3_struc(n)%bgf_fixed     (n_tiles) ) 
       allocate( cmem3_struc(n)%bgf_total     (n_tiles) ) 
       allocate( cmem3_struc(n)%bgf_v     (n_tiles) ) 
       allocate( cmem3_struc(n)%k_lai2vgf     (n_tiles) ) 
       ! snow
       allocate( cmem3_struc(n)%sn_t    (n_tiles) ) 
       allocate( cmem3_struc(n)%sn_moist(n_tiles) ) 
       allocate( cmem3_struc(n)%sn_density(n_tiles) ) 
       allocate( cmem3_struc(n)%sn_depth(n_tiles) ) 
       allocate( cmem3_struc(n)%sn_gsize(n_tiles) ) 
       allocate( cmem3_struc(n)%sn_depth_bias_fac(n_tiles) ) 
       ! output
       allocate( cmem3_struc(n)%em     (n_tiles, n_fs) ) 
       allocate( cmem3_struc(n)%emH     (n_tiles, n_fs) ) 
       allocate( cmem3_struc(n)%emV     (n_tiles, n_fs) ) 
       allocate( cmem3_struc(n)%toaH    (n_tiles, n_fs) ) 
       allocate( cmem3_struc(n)%toaV    (n_tiles, n_fs) ) 

! Initialize values of time invariant parameters/properties
       cmem3_struc(n)%srmax=0.002
       cmem3_struc(n)%srmax2srmin =  0.99   ! % srmin as expressed as fraction of srmax
       cmem3_struc(n)%rho_veg   =  950.0        ! kg/m3, default
       cmem3_struc(n)%m_d       =  0.3          ! -, default
       cmem3_struc(n)%d_leaf    =  0.1          ! mm, default
       cmem3_struc(n)%sn_moist  = 0.1    !  cm3/cm3, default
       cmem3_struc(n)%sn_gsize  = 3.0    !  mm , CMEM default
       cmem3_struc(n)%vwc2lai   = 0.5    !  vwc = 0.5 * lai 
       cmem3_struc(n)%wc_fac    = 1.0    ! 
       cmem3_struc(n)%sn_depth_bias_fac = 1.0  ! multiplicative factor 
       cmem3_struc(n)%k_lai2vgf = 0.52  ! redefined now as exponent for lai-gvf relationship (formerly bare ground fraction) 
       cmem3_struc(n)%bgf_fixed = 0.01  ! default is zero bare ground fraction
    enddo

  end subroutine CMEM3_initialize 


  subroutine add_fields_toState(n, inState,varname)

    use LIS_logMod,   only : LIS_verify
    use LIS_coreMod,  only : LIS_vecTile

    implicit none 

    integer            :: n 
    type(ESMF_State)   :: inState
    character(len=*)   :: varname

    type(ESMF_Field)     :: varField
    type(ESMF_ArraySpec) :: arrspec
    integer              :: status
    real :: sum
    call ESMF_ArraySpecSet(arrspec,rank=1,typekind=ESMF_TYPEKIND_R4,&
         rc=status)
    call LIS_verify(status)

    varField = ESMF_FieldCreate(arrayspec=arrSpec, & 
         grid=LIS_vecTile(n), name=trim(varname), &
         rc=status)
    call LIS_verify(status, 'Error in field_create of '//trim(varname))
    
    call ESMF_StateAdd(inState, (/varField/), rc=status)
    call LIS_verify(status, 'Error in StateAdd of '//trim(varname))

  end subroutine add_fields_toState

  subroutine CMEM3_f2t(n)

    use LIS_coreMod,         only : LIS_rc
    use LIS_metforcingMod,  only : LIS_FORC_State
    use LIS_logMod,          only : LIS_verify
!    USE Profile_Utility, only: sa_to_mr, mr_to_ppmv
!    USE Units_Conversion
    use LIS_FORC_AttributesMod

    implicit none

    integer, intent(in)    :: n 

    integer                :: k, t, kk
    integer                :: status
    type(ESMF_Field)       :: lpressField, pressField, tmpField, q2Field,o3Field
    real,   allocatable        :: lpress(:), press(:), tmp(:), q2(:),o3(:)

  end subroutine CMEM3_f2t


  subroutine CMEM3_geometry(n)
    implicit none
    integer, intent(in)    :: n
  ! do nothing for now
  end subroutine CMEM3_geometry 

  subroutine CMEM3_run(n)
! !USES: 
    use LIS_coreMod, only : LIS_rc, LIS_domain, LIS_surface
    use LIS_RTMMod, only : LIS_sfcState, LIS_forwardState
    use LIS_logMod, only : LIS_verify, LIS_getNextUnitNumber
    use LIS_histDataMod
    
    implicit none

    integer, intent(in) :: n 

    integer             :: status
    integer             :: t, ns, nc, ln, nf, i, j
    real,    pointer    :: wind_speed(:), land_coverage(:), land_temperature(:), snow_temperature(:),&
       soil_moisture_content(:), soil_temperature(:), canopy_water_content(:), &
       vegetation_fraction(:), snow_depth(:), snow_density(:), land_type(:), & 
       snow_coverage(:), leaf_area_index(:), relative_soil_moisture_content(:)
    real, pointer :: tbval(:)

    real :: emV(100), emH(100), toaH, toaV
    real :: bare_emV, bare_emH, bare_toaH, bare_btoaV   ! bare soil 
    real :: veg_emV, veg_emH, veg_toaH, veg_toaV   ! vegeted surface 
    real :: snow_emV, snow_emH, snow_toaH, snow_toaV   ! vegeted surface 
    real :: water_emV, water_emH, water_toaH, water_toaV   ! water bodies 
    real :: sr_temp ! soil moisture-adjusted sr
!    real, parameter :: bgf_fixed=0.5 ! bgf_fixed=bare ground fraction

!   map surface properties to SFC
    call getsfcvar(LIS_sfcState(n), "Wind Speed", wind_speed)
    call getsfcvar(LIS_sfcState(n), "Land Coverage", land_coverage)
    call getsfcvar(LIS_sfcState(n), "Land Type", land_type)
    call getsfcvar(LIS_sfcState(n), "Land Temperature", land_temperature)
    call getsfcvar(LIS_sfcState(n), "Snow Temperature", snow_temperature)
    call getsfcvar(LIS_sfcState(n), "Soil Moisture Content",soil_moisture_content)
    call getsfcvar(LIS_sfcState(n), "Soil Temperature", soil_temperature)
    call getsfcvar(LIS_sfcState(n), "Canopy Water Content", canopy_water_content)
    call getsfcvar(LIS_sfcState(n), "Vegetation Fraction", vegetation_fraction)
    call getsfcvar(LIS_sfcState(n), "Snow Coverage", snow_coverage)
    call getsfcvar(LIS_sfcState(n), "Snow Depth", snow_depth)
    call getsfcvar(LIS_sfcState(n), "Snow Density", snow_density)
    call getsfcvar(LIS_sfcState(n), "Leaf Area Index", leaf_area_index)
    call getsfcvar(LIS_sfcState(n), "Relative Soil Moisture Content", relative_soil_moisture_content)
    
!---------------------------------------------
! Tile loop 
!--------------------------------------------
    do t=1, LIS_rc%ntiles(n)
       cmem3_struc(n)%tsurf(t)     = land_temperature(t)
       cmem3_struc(n)%tau_atm0(t)  = 0.0 
       cmem3_struc(n)%tb_ad(t) 	= 0.0
       cmem3_struc(n)%tb_au(t) 	= 0.0

!  Using soil_temperature results in strong diurnal signal from noah3.2 due to the 
!   considerable difference
!       cmem3_struc(n)%tsoil(t)     = soil_temperature(t)
       cmem3_struc(n)%tsoil(t)     = land_temperature(t)

       cmem3_struc(n)%wc(t)  = soil_moisture_content(t) 
! Soil moisture correction: rescaling method
       i = LIS_domain(n)%tile(t)%col
       j = LIS_domain(n)%tile(t)%row
       If (sm_correction_struc(n)%c_type .eq. 1) then
          cmem3_struc(n)%wc(t)  =  & 
               sm_correction_struc(n)%dst_mean(i, j) + &
               (soil_moisture_content(t)-&
               sm_correction_struc(n)%src_mean(i, j) ) * &
               sm_correction_struc(n)%dst_sigma(i, j)/&
               sm_correction_struc(n)%src_sigma(i, j)
          
          if ( cmem3_struc(n)%wc(t) .LT. 0 ) &
               cmem3_struc(n)%wc(t) = 0.0 
       End if
       cmem3_struc(n)%wc(t)  = cmem3_struc(n)%wc(t) * &
            cmem3_struc(n)%wc_fac(t) ! m3/m3
       cmem3_struc(n)%relsmc(t)= relative_soil_moisture_content(t)       
       !      cmem3_struc(n)%sr(t) 	= 0.002                       ! m 
       cmem3_struc(n)%sr(t) = cmem3_struc(n)%srmax(t)  &
            *((cmem3_struc(n)%srmax2srmin(t)-1) &
            * cmem3_struc(n)%relsmc(t) +1)

       cmem3_struc(n)%sand(t) 	= LIS_surface(n,&
            LIS_rc%lsm_index)%tile(t)%sand*100.0
       cmem3_struc(n)%clay(t) 	= LIS_surface(n,&
            LIS_rc%lsm_index)%tile(t)%clay*100.0
       cmem3_struc(n)%vtype(t)     = nint(land_type(t))   ! asssuming UMD 
       cmem3_struc(n)%tveg(t)      = land_temperature(t)
       cmem3_struc(n)%h_veg(t)     =  10.0           ! not used yet
!       cmem3_struc(n)%veg_frac(t)  = vegetation_fraction(t)
       IF (leaf_area_index(t) .eq. 0.00) THEN
           leaf_area_index(t) = 0.01
       END IF
       cmem3_struc(n)%bgf_v(t)      =  1.0*exp(-1.0*cmem3_struc(n)%k_lai2vgf(t)*(leaf_area_index(t)/(1.0-cmem3_struc(n)%bgf_fixed(t))))
       cmem3_struc(n)%bgf_total(t)  =  cmem3_struc(n)%bgf_fixed(t) + cmem3_struc(n)%bgf_v(t)*(1.0-cmem3_struc(n)%bgf_fixed(t))
       cmem3_struc(n)%lai(t)       =  leaf_area_index(t)/(1.0-cmem3_struc(n)%bgf_total(t))

!if (t.eq.1) then
!   print *, cmem3_struc(n)%bgf_total(t), cmem3_struc(n)%bgf_fixed(t), cmem3_struc(n)%bgf_v(t), &
!        cmem3_struc(n)%k_lai2vgf(t), leaf_area_index(t),  cmem3_struc(n)%lai(t)
!endif

!       cmem3_struc(n)%lai(t)       =  leaf_area_index(t)/cmem3_struc(n)%veg_frac(t)
!      cmem3_struc(n)%vwc2lai      =  0.5                   ! used to est vwc 
!      cmem3_struc(n)%rho_veg(t)   =  950.0        ! kg/m3, default
!      cmem3_struc(n)%m_d(t)       =  0.3          ! -, default
!      cmem3_struc(n)%d_leaf(t)    =  0.1          ! mm, default

       cmem3_struc(n)%sn_t(t)	= snow_temperature(t)
!      cmem3_struc(n)%sn_moist(t)  = 0.1    !  cm3/cm3, default 
       cmem3_struc(n)%sn_density(t)= snow_density(t)    ! g/cm3
       cmem3_struc(n)%sn_depth(t)  = snow_depth(t) * &
            cmem3_struc(n)%sn_depth_bias_fac(t)      ! m

!      cmem3_struc(n)%sn_gsize(t)  = 3.0    !  mm , CMEM default

! -----------------------------------------------------------------------------
! New formulation: canopy water content = canopy interception (cmc) + vegetation
!  water content (vwc) ( and + moisture in litter, possibly in future).
! VWC here is estimated from LAI, with the simple formulation VWC (mm) = 0.5 * LAI
! which has been used in CRTM and L-MEB, but field experiments show that such
! a relation depends on vegetation type (e.g., Anderson et al., 2004, Remote Sens. 
! Envi. 92, 447-464.)
! -----------------------------------------------------------------------------

       cmem3_struc(n)%wc_veg(t)    =  canopy_water_content(t) + & 
            leaf_area_index(t) * cmem3_struc(n)%vwc2lai(t) 

       emH =0.0
       emV =0.0
       ! frequency loop 
       Do nf = 1, cmem3_struc(n)%nfs
          bare_emH=0.0
          bare_emV=0.0
          
          veg_emH=0.0
          veg_emV=0.0
	
          snow_emH=0.0
          snow_emV=0.0
          
          water_emH=0.0
          water_emV=0.0

! YDT: 4/10/2012. CMEM's snow module tends to crash. Disable it for now
! by resetting snow_coverage to 0, and land_coverage to 1. 
          snow_coverage(t) = 0.0
          land_coverage(t) = 1.0

       ! land_coverage + snow_coverage = 1
          if (snow_coverage(t) .GT. 0) then 
             cmem3_struc(n)%vtype(t)     = nint(land_type(t))   ! asssuming UMD 
             call CMEM3_Compute(n, t, nf) 
             snow_emH = cmem3_struc(n)%emH(t, nf) 
             snow_emV = cmem3_struc(n)%emV(t, nf) 
!         write(*, '(A, 4F10.3)')"snow_depth:fghz:emH:emV ", snow_depth(t), cmem3_struc(n)%fghz(nf),  &
!                    snow_emH, snow_emV 
          end if
          if (land_coverage(t) .GT. 0) then      ! snow-free land 
             cmem3_struc(n)%sn_depth(t)   = 0.0   
             
!             if (cmem3_struc(n)%veg_frac(t) .GT. 0 ) then  ! calculate veg'd surface
             if ((1.0-cmem3_struc(n)%bgf_total(t)) .GT. 0 ) then  ! calculate veg'd surface
                cmem3_struc(n)%vtype(t)     = nint(land_type(t))   ! asssuming UMD 
                call CMEM3_Compute(n, t, nf) 
                veg_emH = cmem3_struc(n)%emH(t, nf) 
                veg_emV = cmem3_struc(n)%emV(t, nf) 
             end if
             
             cmem3_struc(n)%vtype(t) = 12   ! bare ground for UMD
             call CMEM3_Compute(n, t, nf) 
             bare_emH = cmem3_struc(n)%emH(t, nf) 
             bare_emV = cmem3_struc(n)%emV(t, nf) 
             
          end if
       ! radiation tile averaging   TEMP TESTING OF IMPACT OF BARE GROUND 
          cmem3_struc(n)%emH(t, nf) = snow_emH * snow_coverage(t) +  &
               veg_emH * (1.0-cmem3_struc(n)%bgf_total(t)) * land_coverage(t) + &
               bare_emH * (cmem3_struc(n)%bgf_total(t)) * land_coverage(t)
          cmem3_struc(n)%emV(t, nf) = snow_emV * snow_coverage(t) +  &
               veg_emV * (1.0-cmem3_struc(n)%bgf_total(t)) * land_coverage(t) + &
               bare_emV * (cmem3_struc(n)%bgf_total(t)) * land_coverage(t) 

!!!          cmem3_struc(n)%emH(t, nf) = snow_emH * snow_coverage(t) +  &
!!!               veg_emH * (1.0-cmem3_struc(n)%bgf_fixed(t)-cmem3_struc(n)%veg_frac(t)) * land_coverage(t) + &
!!!               bare_emH * (cmem3_struc(n)%bgf_fixed(t)+cmem3_struc(n)%veg_frac(t)) * land_coverage(t)
!!!          cmem3_struc(n)%emV(t, nf) = snow_emV * snow_coverage(t) +  &
!!!               veg_emV * (1.0-cmem3_struc(n)%bgf_fixed(t)) * land_coverage(t) + &
!!!               bare_emV * (cmem3_struc(n)%bgf_fixed(t) ) * land_coverage(t) 
          
!!!!          cmem3_struc(n)%emH(t, nf) = snow_emH * snow_coverage(t) +  &
!!!!               veg_emH * (0.0) * land_coverage(t) + &
!!!!               bare_emH * (1.0) * land_coverage(t)
!!!!          cmem3_struc(n)%emV(t, nf) = snow_emV * snow_coverage(t) +  &
!!!!               veg_emV * (0.0) * land_coverage(t) + &
!!!!               bare_emV * (1.0 ) * land_coverage(t) 
          
          ! veg tile averaging, skipping for now
       !write(*, *)"nf=", nf, " emH=", emH(nf), " emV=", emV(nf) 
       !emH(nf) = emH(nf) + cmem3_struc(n)%emH(t, nf) * LIS_domain(n)%tile(t)%fgrd
       !emV(nf) = emV(nf) + cmem3_struc(n)%emV(t, nf) * LIS_domain(n)%tile(t)%fgrd
       
          call LIS_diagnoseRTMOutputVar(n, t,LIS_MOC_RTM_EMISSIVITY,value=    &
               cmem3_struc(n)%emH(t,nf),vlevel=(nf-1)*2+1,unit="-",&
               direction="-")
          call LIS_diagnoseRTMOutputVar(n, t,LIS_MOC_RTM_EMISSIVITY,value=    &
               cmem3_struc(n)%emV(t,nf),vlevel=(nf-1)*2+2,unit="-",&
               direction="-")

       
#if 0 
          call LIS_diagnoseRTMOutputVar(n, t,LIS_MOC_RTM_TB,value=    &
               cmem3_struc(n)%emH(t,nf)*land_temperature(t),&
               vlevel=(nf-1)*2+1,unit="K",direction="-")
          call LIS_diagnoseRTMOutputVar(n, t,LIS_MOC_RTM_TB,value=    &
               cmem3_struc(n)%emV(t,nf)*land_temperature(t),&
               vlevel=(nf-1)*2+2,unit="K",direction="-")
#endif

          call LIS_diagnoseRTMOutputVar(n, t,LIS_MOC_RTM_TB,value=    &
               cmem3_struc(n)%toaH(t,nf),vlevel=(nf-1)*2+1,unit="K",&
               direction="-")
          call LIS_diagnoseRTMOutputVar(n, t,LIS_MOC_RTM_TB,value=    &
               cmem3_struc(n)%toaV(t,nf),vlevel=(nf-1)*2+2,unit="K",&
               direction="-")
          
       end do ! frequency loop
       

       call LIS_diagnoseRTMOutputVar(n, t, LIS_MOC_RTM_SM, value=       &
            real( cmem3_struc(n)%wc(t) ),             &
            vlevel=1, unit="m3/m3",direction="-")
       
    End Do  ! tile loop

    do j=1, cmem3_struc(n)%nfs
       call getsfcvar(LIS_forwardState(n),&
            trim(cmem3_struc(n)%chName(j))//"_Hpol", tbval)
       tbval = cmem3_struc(n)%toaH(:,j)

       call getsfcvar(LIS_forwardState(n),&
            trim(cmem3_struc(n)%chName(j))//"_Vpol", tbval)
       tbval = cmem3_struc(n)%toaV(:,j)
    enddo
  end subroutine CMEM3_run
  

  subroutine CMEM3_Compute ( n, t, nf )

! !USES:
    use LIS_coreMod, only : LIS_rc
    use LIS_RTMMod, only : LIS_sfcState
    use LIS_logMod, only : LIS_verify
    implicit none

    integer, intent(in) :: n, t, nf
    integer :: i

! Debugging output
!if (t .eq. 1) then 
if (t .eq. 0) then 
   write(*, *)"tscount=", LIS_rc%tscount 
   write(*, '(A40, F12.4)')"cmem3_struc(n)%fghz(nf)", cmem3_struc(n)%fghz(nf)
   write(*, '(A40, F12.4)')"cmem3_struc(n)%theta(nf)", cmem3_struc(n)%theta(nf)
   write(*, '(A40, F12.4)')"cmem3_struc(n)%tsurf(t)", cmem3_struc(n)%tsurf(t)
   write(*, '(A40, F12.4)')"cmem3_struc(n)%tsoil(t)", cmem3_struc(n)%tsoil(t)
   write(*, '(A40, F12.4)')"cmem3_struc(n)%wc(t)", cmem3_struc(n)%wc(t)
   write(*, '(A40, F12.4)')"cmem3_struc(n)%sr(t)", cmem3_struc(n)%sr(t)
   write(*, '(A40, F12.4)')"cmem3_struc(n)%srmax2srmin(t)", cmem3_struc(n)%srmax2srmin(t)
   write(*, '(A40, F12.4)')"cmem3_struc(n)%relsmc(t)", cmem3_struc(n)%relsmc(t)
   write(*, '(A40, F12.4)')"cmem3_struc(n)%sand(t)", cmem3_struc(n)%sand(t)
   write(*, '(A40, F12.4)')"cmem3_struc(n)%clay(t)", cmem3_struc(n)%clay(t)
   write(*, '(A40, F12.4)')"cmem3_struc(n)%vtype(t)", cmem3_struc(n)%vtype(t)
   write(*, '(A40, F12.4)')"cmem3_struc(n)%tveg(t)", cmem3_struc(n)%tveg(t)
   write(*, '(A40, F12.4)')"cmem3_struc(n)%h_veg(t)", cmem3_struc(n)%h_veg(t)
   write(*, '(A40, F12.4)')"cmem3_struc(n)%lai(t)", cmem3_struc(n)%lai(t)
   write(*, '(A40, F12.4)')"cmem3_struc(n)%wc_veg(t)", cmem3_struc(n)%wc_veg(t)
   write(*, '(A40, F12.4)')"cmem3_struc(n)%rho_veg(t)", cmem3_struc(n)%rho_veg(t)
   write(*, '(A40, F12.4)')"cmem3_struc(n)%m_d(t)", cmem3_struc(n)%m_d(t)
   write(*, '(A40, F12.4)')"cmem3_struc(n)%d_leaf(t)", cmem3_struc(n)%d_leaf(t)
   write(*, '(A40, F12.4)')"cmem3_struc(n)%sn_t(t)", cmem3_struc(n)%sn_t(t)
   write(*, '(A40, F12.4)')"cmem3_struc(n)%sn_moist(t)", cmem3_struc(n)%sn_moist(t)
   write(*, '(A40, F12.4)')"cmem3_struc(n)%sn_density(t)", cmem3_struc(n)%sn_density(t)
   write(*, '(A40, F12.4)')"cmem3_struc(n)%sn_depth(t)", cmem3_struc(n)%sn_depth(t)
   write(*, '(A40, F12.4)')"cmem3_struc(n)%sn_gsize(t)", cmem3_struc(n)%sn_gsize(t)
endif

    call lis_mem ( cmem3_struc(n)%fghz(nf),		& 
                   cmem3_struc(n)%theta(nf),	& 
                   cmem3_struc(n)%tsurf(t),	&  ! atmos
                   cmem3_struc(n)%tau_atm0(t),	&
                   cmem3_struc(n)%tb_ad(t),	& 
                   cmem3_struc(n)%tb_au(t),	& 
                   cmem3_struc(n)%tsoil(t),	&  ! soil
                   cmem3_struc(n)%wc(t),   	& 
                   cmem3_struc(n)%sr(t),   	& 
                   cmem3_struc(n)%sand(t), 	& 
                   cmem3_struc(n)%clay(t), 	& 
                   cmem3_struc(n)%vtype(t),	&  ! veg
                   cmem3_struc(n)%tveg(t),	& 
                   cmem3_struc(n)%h_veg(t),	& 
                   cmem3_struc(n)%lai(t),	& 
                   cmem3_struc(n)%wc_veg(t),	&
                   cmem3_struc(n)%rho_veg(t),	&
                   cmem3_struc(n)%m_d(t),    	&
                   cmem3_struc(n)%d_leaf(t), 	&
                   cmem3_struc(n)%sn_t(t),   	&  ! snow
                   cmem3_struc(n)%sn_moist(t),   	&  
                   cmem3_struc(n)%sn_density(t),  	&  
                   cmem3_struc(n)%sn_depth(t),    	&  
                   cmem3_struc(n)%sn_gsize(t),    	&  
                   cmem3_struc(n)%emH(t, nf),    	&  ! output
                   cmem3_struc(n)%emV(t, nf),    	&  
                   cmem3_struc(n)%toaH(t, nf),   	&  
                   cmem3_struc(n)%toaV(t, nf)      ) 

  end subroutine CMEM3_Compute

  subroutine getsfcvar(sfcState, varname, var)
! !USES: 
    use LIS_logMod,  only : LIS_verify
    
    implicit none
    
    type(ESMF_State)      :: sfcState
    type(ESMF_Field)      :: varField
    character(len=*)      :: varname
    real, pointer         :: var(:)
    integer               :: status

    call ESMF_StateGet(sfcState, trim(varname), varField, rc=status)
    call LIS_verify(status, 'Error in StateGet: CMEM3_handlerMod '//trim(varname))
    call ESMF_FieldGet(varField, localDE=0,farrayPtr=var, rc=status)
    call LIS_verify(status, 'Error in FieldGet: CMEM3_handlerMod '//trim(varname))

  end subroutine getsfcvar


!!$  integer function CMEM3_landmatch(i, classification)
!!$
!!$  !USES:      
!!$    use CRTM_Module
!!$
!!$    implicit none
!!$  ! !ARGUMENTS: 
!!$    integer, intent(in) :: i               ! index into CMEM3 lookup table
!!$    integer, intent(in) :: classification  ! classification scheme, e.g., UMD
!!$  !                                          (1), USGS (2), UMD (3), IGBP (4)
!!$  ! Description
!!$  ! Output 
!!$    integer :: j                           ! temp variable to store returned value
!!$
!!$  ! UMD Land Cover Classification from:  http://www.geog.umd.edu/landcover/1km-map/meta-data.html
!!$  !kwh: ARBITRARY VALUES CURRENTLY ASSIGNED FOR TESTING.
!!$
!!$    integer, dimension (0:13), parameter :: umd_cmem3_match = &
!!$      (/  &
!!$           INVALID_LAND              , &        ! 0  - WATER (and Goode's interrupted space)   
!!$           PINE_FOREST               , &        ! 1  - EVERGREEN NEEDLELEAF FOREST             
!!$           BROADLEAF_FOREST          , &        ! 2  - EVERGREEN BROADLEAF FOREST              
!!$           PINE_FOREST               , &        ! 3  - DECIDUOUS NEEDLELEAF FOREST             
!!$           BROADLEAF_FOREST          , &        ! 4  - DECIDUOUS BROADLEAF FOREST              
!!$           BROADLEAF_PINE_FOREST     , &        ! 5  - MIXED FOREST                            
!!$           BROADLEAF_BRUSH           , &        ! 6  - WOODLAND                                
!!$           BROADLEAF_BRUSH           , &        ! 7  - WOODED GRASSLAND                        
!!$           SCRUB                     , &        ! 8  - CLOSED SHRUBLAND                        
!!$           SCRUB_SOIL                , &        ! 9  - OPEN SHRUBLAND                          
!!$           GRASS_SCRUB               , &        ! 10 - GRASSLAND                               
!!$           TILLED_SOIL               , &        ! 11 - CROPLAND                                
!!$           COMPACTED_SOIL            , &        ! 12 - BARE GROUND                             
!!$           URBAN_CONCRETE            &          ! 13 - URBAN AND BUILT-UP                      
!!$       /)                                                         
!!$               
!!$  !                
!!$  !USGS Land Use/Land Cover System Legend (Modified Level 2) (as taken from USGSEDC)
!!$  !kwh: NEED TO PUT IN CMEM3 MATCH HERE, AS WITH GFS.  CURRENLT ARBITRARY VALUES ASSIGNED
!!$
!!$    integer, dimension (24), parameter :: usgs_cmem3_match = &
!!$    (/  &                                         !Value 	Code 	Description
!!$       URBAN_CONCRETE            , &          	! 1  	100 	Urban and Built-Up Land
!!$       TILLED_SOIL               , &          	! 2  	211 	Dryland Cropland and Pasture
!!$       TILLED_SOIL               , &          	! 3  	212 	Irrigated Cropland and Pasture
!!$       TILLED_SOIL               , &          	! 4  	213 	Mixed Dryland/Irrigated Cropland and Pasture
!!$       TILLED_SOIL               , &          	! 5  	280 	Cropland/Grassland Mosaic
!!$       TILLED_SOIL               , &          	! 6  	290 	Cropland/Woodland Mosaic
!!$       MEADOW_GRASS              , &          	! 7  	311 	Grassland
!!$       SCRUB                     , &          	! 8  	321 	Shrubland
!!$       GRASS_SCRUB               , &          	! 9  	330 	Mixed Shrubland/Grassland
!!$       BROADLEAF_BRUSH           , &          	! 10 	332 	Savanna
!!$       BROADLEAF_FOREST          , &          	! 11 	411 	Deciduous Broadleaf Forest
!!$       PINE_FOREST               , &          	! 12 	412 	Deciduous Needleleaf Forest
!!$       BROADLEAF_FOREST          , &          	! 13 	421 	Evergreen Broadleaf Forest
!!$       PINE_FOREST               , &         	! 14 	422 	Evergreen Needleleaf Forest
!!$       BROADLEAF_PINE_FOREST     , &         	! 15 	430 	Mixed Forest
!!$       INVALID_LAND              , &		! 16 	500 	Water Bodies
!!$       BROADLEAF_BRUSH           , &		! 17 	620 	Herbaceous Wetland
!!$       BROADLEAF_BRUSH           , &		! 18 	610 	Wooded Wetland
!!$       COMPACTED_SOIL 	       , &		! 19 	770 	Barren or Sparsely Vegetated
!!$       TUNDRA		       , &		! 20 	820 	Herbaceous Tundra
!!$       TUNDRA  		       , &		! 21 	810 	Wooded Tundra
!!$       TUNDRA  		       , &		! 22 	850 	Mixed Tundra
!!$       TUNDRA  		       , &		! 23 	830 	Bare Ground Tundra
!!$       INVALID_LAND   	       &		! 24 	900 	Snow or Ice  
!!$     /)  
!!$                                  
!!$  !  integer, parameter :: USGS_  =99 	  	! Interrupted Areas (Goodes Homolosine Projection)
!!$  !  integer, parameter :: USGS_  =100 	Missing Data
!!$  !
!!$  !   GFS/GDAS vegtype taken from Weizhong Zeng driver
!!$
!!$    integer, dimension (13), parameter :: gfs_cmem3_match = &
!!$    (/  &                              
!!$         BROADLEAF_FOREST         , &   ! 1  - BROADLEAF-EVERGREEN TREES  (TROPICAL FOREST)
!!$         BROADLEAF_FOREST         , &   ! 2  - BROADLEAF-DECIDUOUS TREES
!!$         BROADLEAF_PINE_FOREST    , &   ! 3  - BROADLEAF AND NEEDLELEAF TREES (MIXED FOREST)
!!$         PINE_FOREST              , &   ! 4  - NEEDLELEAF-EVERGREEN TREES
!!$         PINE_FOREST              , &   ! 5  - NEEDLELEAF-DECIDUOUS TREES (LARCH)
!!$         BROADLEAF_BRUSH          , &   ! 6  - BROADLEAF TREES WITH GROUNDCOVER (SAVANNA)
!!$         SCRUB                    , &   ! 7  - GROUNDCOVER ONLY (PERENNIAL)  
!!$         SCRUB                    , &   ! 8  - BROADLEAF SHRUBS WITH PERENNIAL GROUNDCOVER
!!$         SCRUB_SOIL               , &   ! 9  - BROADLEAF SHRUBS WITH BARE SOIL
!!$         TUNDRA                   , &   ! 10 - DWARF TREES AND SHRUBS WITH GROUNDCOVER (TUNDRA)
!!$         COMPACTED_SOIL           , &   ! 11 - BARE SOIL
!!$         TILLED_SOIL              , &   ! 12 - CULTIVATIONS (THE SAME PARAMETERS AS FOR TYPE 7)
!!$         COMPACTED_SOIL           &     ! 13 - GLACIAL (THE SAME PARAMETERS AS FOR TYPE 11) 
!!$    /)                                
!!$
!!$
!!$  ! IGBP Land Cover Categories as taken from USGS EDC
!!$  !kwh: NEED TO PUT IN CMEM3 MATCH HERE, AS WITH GFS.  CURRENTLY ARBITRARY VALUES ASSIGNED
!!$    integer, dimension (17), parameter :: igbp_cmem3_match = &
!!$    (/  &                              
!!$         PINE_FOREST                , &               ! 1  Evergreen Needleleaf Forest
!!$         BROADLEAF_FOREST           , &               ! 2  Evergreen Broadleaf Forest
!!$         PINE_FOREST                , &               ! 3  Deciduous Needleleaf Forest
!!$         BROADLEAF_FOREST           , &               ! 4  Deciduous Broadleaf Forest
!!$         BROADLEAF_PINE_FOREST      , &               ! 5  Mixed Forest
!!$         SCRUB                      , &               ! 6  Closed Shrublands
!!$         SCRUB_SOIL                 , &               ! 7  Open Shrublands
!!$         BROADLEAF_BRUSH            , &               ! 8  Woody Savannas
!!$         BROADLEAF_BRUSH            , &               ! 9  Savannas
!!$         MEADOW_GRASS               , &               ! 10 Grasslands
!!$         BROADLEAF_BRUSH            , &               ! 11 Permanent Wetlands
!!$         TILLED_SOIL                , &               ! 12 Croplands
!!$         URBAN_CONCRETE             , &               ! 13 Urban and Built-Up
!!$         TILLED_SOIL                , &               ! 14 Cropland/Natural Vegetation Mosaic
!!$         INVALID_LAND               , &               ! 15 Snow and Ice
!!$         COMPACTED_SOIL             , &               ! 16 Barren or Sparsely Vegetated
!!$         INVALID_LAND               &                 ! 17 Water Bodies                                     
!!$  !                                                   ! 99 Interrupted Areas (Goodes Homolosine Projection)
!!$  !                                                   !100 Missing Data           
!!$    /)                                
!!$                                    
!!$  ! 
!!$  ! !DESCRIPTION
!!$  ! scheme:
!!$  ! #1 use the UMD landcover
!!$  ! #2 use the USGS landcover data
!!$  ! #3 use the GFS landcover data
!!$  ! #4 use the IGBP landcover data
!!$  !
!!$  ! This subroutine matches the land (veg type) classification
!!$  ! to that of CMEM3
!!$
!!$  !EOP
!!$
!!$  if     (classification .eq. 1) then !UMD
!!$    !need to bounce matchings off someone; similar matchings to weizhong driver
!!$    j = umd_cmem3_match(i)
!!$  elseif (classification .eq. 2) then !USGS
!!$    !insert error code for 'not yet implemented'
!!$    !j = usgs_cmem3_match(i)
!!$  elseif (classification .eq. 3) then !GFS
!!$    j = gfs_cmem3_match(i)
!!$  elseif (classification .eq. 4) then !IGBP
!!$    !insert error code for 'not yet implemented'
!!$    !j = igbp_cmem3_match(i)
!!$  else
!!$    !kwh: insert error code
!!$  end if
!!$  CMEM3_landmatch = j
!!$  return
!!$end function  CMEM3_landmatch

#endif
end module CMEM3_Mod



