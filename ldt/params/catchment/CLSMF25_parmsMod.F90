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
module CLSMF25_parmsMod
!BOP
!
! !MODULE: CLSMF25_parmsMod
!
! !DESCRIPTION:
!  The code in this file implements routines to read Catchment 
!  parameter datasets.
!
!  \subsubsection{Overview}
!  This routines in this module provides routines to read the 
!  Catchment LSM  data and allows the users to 
!  specify the frequency of climatology (in months). 
!  The climatological data is temporally interpolated  
!  between months to the current simulation date. 
!
! !REVISION HISTORY:
!
!  08 Aug 2005: Sujay Kumar; Initial implementation
!  20 Nov 2012: K. Arsenault: Updated for Catchment LSM parameters.
!
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif
  use ESMF
  use LDT_coreMod
  use LDT_historyMod
  use LDT_gfracMod
  use LDT_albedoMod
  use LDT_paramDataMod
  use LDT_logMod
  use LDT_constantsMod, only : LDT_CONST_PATH_LEN

  implicit none

  PRIVATE

!------------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!------------------------------------------------------------------------------
  public :: catchmentParms_init    !allocates memory for required structures
  public :: catchmentParms_writeHeader
  public :: catchmentParms_writeData

!------------------------------------------------------------------------------
! !PUBLIC TYPES:
!------------------------------------------------------------------------------
  public :: CLSMF25_struc

  type, public :: clsmf25_type_dec

     character(len=LDT_CONST_PATH_LEN) :: albnirfile
     character(len=LDT_CONST_PATH_LEN) :: albvisfile
     real          :: catchparms_gridDesc(20)
     character*50  :: catchparms_proj
     character*50  :: catchparms_gridtransform
!  - Catchment (F2.5):
!     character*100 :: modisdir        ! sub-directory for MODIS files
!     character*140 :: tile_coord_file ! tile coordinate file
!     character*140 :: tile_veg_file   ! tile vegetation file
     character(len=LDT_CONST_PATH_LEN) :: soilparamfile   ! soil parameters file
     character(len=LDT_CONST_PATH_LEN) :: sltsfile        ! surface layer timescales file
     character(len=LDT_CONST_PATH_LEN) :: topo_ar_file    ! topography parameters file
     character(len=LDT_CONST_PATH_LEN) :: topo_bf_file    ! topography parameters file
     character(len=LDT_CONST_PATH_LEN) :: topo_ts_file    ! topography parameters file
     character(len=LDT_CONST_PATH_LEN) :: catchgreenfile  ! greenness climatology file
     character(len=LDT_CONST_PATH_LEN) :: catchlaifile    ! LAI climatology file
     real          :: dzsfcrd         ! CLSM top soil layer depth from ldt.config
     real          :: addbdrckcrd     ! CLSM add to bedrock depth from ldt.config

     type(LDT_paramEntry) :: psisat      ! saturated soil moisture potential
     type(LDT_paramEntry) :: bexp        ! Clapp-Hornberger parameter
     type(LDT_paramEntry) :: wpwet       ! wilting point wetness
     type(LDT_paramEntry) :: bdrckdpth   ! depth to bedrock 
     type(LDT_paramEntry) :: ksat        ! saturated hydraulic conductivity (SATDK; m s-1)
     type(LDT_paramEntry) :: gnu         ! vertical decay factor for transmissivity
     type(LDT_paramEntry) :: ars1        ! Wetness parameters
     type(LDT_paramEntry) :: ars2         
     type(LDT_paramEntry) :: ars3
     type(LDT_paramEntry) :: ara1        ! Shape parameters
     type(LDT_paramEntry) :: ara2
     type(LDT_paramEntry) :: ara3
     type(LDT_paramEntry) :: ara4
     type(LDT_paramEntry) :: arw1        ! Minimum Theta parameters:
     type(LDT_paramEntry) :: arw2
     type(LDT_paramEntry) :: arw3
     type(LDT_paramEntry) :: arw4
     type(LDT_paramEntry) :: bf1         ! Baseflow topographic params
     type(LDT_paramEntry) :: bf2
     type(LDT_paramEntry) :: bf3
     type(LDT_paramEntry) :: tsa1        ! Water transfer parameters
     type(LDT_paramEntry) :: tsa2
     type(LDT_paramEntry) :: tsb1        
     type(LDT_paramEntry) :: tsb2
     type(LDT_paramEntry) :: atau
     type(LDT_paramEntry) :: btau
     type(LDT_paramEntry) :: albnirdir   ! Albedo NIR direct scale factor (CLSM F2.5)
     type(LDT_paramEntry) :: albnirdif   ! Albedo NIR diffuse scale factor (CLSM F2.5)
     type(LDT_paramEntry) :: albvisdir   ! Albedo VIS direct scale factor (CLSM F2.5)
     type(LDT_paramEntry) :: albvisdif   ! Albedo VIS diffuse scale factor (CLSM F2.5)

  end type clsmf25_type_dec

  type(clsmf25_type_dec), allocatable :: CLSMF25_struc(:)

contains

!BOP
! 
! !ROUTINE: catchmentParms_init
! \label{catchmentParms_init}
! 
! !INTERFACE:
  subroutine catchmentParms_init(flag)
! !USES:
    use LDT_fileIOMod, only : LDT_readDomainConfigSpecs
    use LDT_logMod,    only : LDT_verify
    use LDT_paramOptCheckMod, only: &! LDT_catchparmsOptChecks, &
                       LDT_gridOptChecks
! !DESCRIPTION:
!
! Allocates memory for data structures for reading 
! the Catchment LSM Parameter datasets.
! 
!  The routines invoked are: 
!  \begin{description}
!   \item[catchmentParms](\ref{catchmentParms}) \newline
!    calls the registry to invoke the Catchment parms setup methods. 
!  \end{description}
!
!EOP
    implicit none
    integer      :: flag
    integer      :: n
    integer      :: rc
    character*50 :: catchparms_proj
    real         :: catchparms_gridDesc(LDT_rc%nnest, 20)
! ________________________________________________________

    write(LDT_logunit,*)" - - - - - - - - - Catchment LSM Parameters - - - - - - - - - - - -"
    
    allocate(CLSMF25_struc(LDT_rc%nnest))

    do n=1,LDT_rc%nnest

       call set_param_attribs(CLSMF25_struc(n)%bexp,"BEXP", &
            full_name="CLSMF2.5 Bexp Clapp-Hornberger parameter")

       call set_param_attribs(CLSMF25_struc(n)%psisat,"PSISAT", &
            full_name="CLSMF2.5 saturated soil moisture potential")

       call set_param_attribs(CLSMF25_struc(n)%wpwet,"WPWET", &
            full_name="CLSMF2.5 wilting point wetness")

       call set_param_attribs(CLSMF25_struc(n)%ksat,"KSAT",&
            units="ms-1",full_name="CLSMF2.5 saturated hydraulic conductivity")

       call set_param_attribs(CLSMF25_struc(n)%gnu,"GNUCLSM",&
            units="m-1", &
            full_name="CLSMF2.5 vertical transm. decay term")

       call set_param_attribs(CLSMF25_struc(n)%ars1,"ARS1CLSM",&
            units="m2kg-1", &
            full_name="CLSMF2.5 (ARS) wetness parameters" )

       call set_param_attribs(CLSMF25_struc(n)%ara1,"ARA1CLSM", &
            units="m2kg-1", &
            full_name="CLSMF2.5 (ARA) topographic shape parameters" )

       call set_param_attribs(CLSMF25_struc(n)%arw1,"ARW1CLSM",&
            units="m2kg-1", &
            full_name="CLSMF2.5 (ARW) minimum theta parameters" )

       call set_param_attribs(CLSMF25_struc(n)%bf1,"BF1CLSM",&
            units="kgm-4", &
            full_name="CLSMF2.5 (BF) baseflow topographic parameters" )

       call set_param_attribs(CLSMF25_struc(n)%tsa1,"TSA1CLSM",&
            full_name="CLSMF2.5 (TS) water transfer parameters" )

       call set_param_attribs(CLSMF25_struc(n)%atau,"ATAUCLSM",&
            full_name="CLSMF2.5 (TAU) topographic tau parameters" )

       call set_param_attribs(CLSMF25_struc(n)%bdrckdpth,"BEDROCKDEPTH",&
            units="mm", &
            full_name="CLSMF2.5 depth to bedrock" )

       call set_param_attribs(CLSMF25_struc(n)%albnirdir,"ALBNIRDIR",&
            vlevels=12, &
            full_name="CLSMF2.5 alb near-IR (direct) scale factor" )

       call set_param_attribs(CLSMF25_struc(n)%albvisdir,"ALBVISDIR",&
            vlevels=12, &
            full_name="CLSMF2.5 alb visible (direct) scale factor" )

       call set_param_attribs(CLSMF25_struc(n)%albnirdif,"ALBNIRDIF",&
            vlevels=12, &
            full_name="CLSMF2.5 alb near-IR (diffuse) scale factor" )

       call set_param_attribs(CLSMF25_struc(n)%albvisdif,"ALBVISDIF",&
            vlevels=12, &
            full_name="CLSMF2.5 alb near-IR (diffuse) scale factor" )

       allocate(CLSMF25_struc(n)%ksat%value(&
            LDT_rc%lnc(n),LDT_rc%lnr(n),&
            CLSMF25_struc(n)%ksat%num_bins))
       
       allocate(CLSMF25_struc(n)%psisat%value(&
            LDT_rc%lnc(n),LDT_rc%lnr(n),&
            CLSMF25_struc(n)%psisat%num_bins))
       
       allocate(CLSMF25_struc(n)%bexp%value(&
            LDT_rc%lnc(n),LDT_rc%lnr(n),&
            CLSMF25_struc(n)%bexp%num_bins))
       
       allocate(CLSMF25_struc(n)%wpwet%value(&
            LDT_rc%lnc(n),LDT_rc%lnr(n),&
            CLSMF25_struc(n)%wpwet%num_bins))
       allocate(CLSMF25_struc(n)%bdrckdpth%value(&
            LDT_rc%lnc(n),LDT_rc%lnr(n),&
            CLSMF25_struc(n)%bdrckdpth%num_bins))
       allocate(CLSMF25_struc(n)%gnu%value(&
            LDT_rc%lnc(n),LDT_rc%lnr(n),&
            CLSMF25_struc(n)%gnu%vlevels))

       allocate(CLSMF25_struc(n)%albnirdir%value(&
            LDT_rc%lnc(n),LDT_rc%lnr(n),&
            CLSMF25_struc(n)%albnirdir%vlevels))
       allocate(CLSMF25_struc(n)%albvisdir%value(&
            LDT_rc%lnc(n),LDT_rc%lnr(n),&
            CLSMF25_struc(n)%albvisdir%vlevels))
       allocate(CLSMF25_struc(n)%albnirdif%value(&
            LDT_rc%lnc(n),LDT_rc%lnr(n),&
            CLSMF25_struc(n)%albnirdif%vlevels))
       allocate(CLSMF25_struc(n)%albvisdif%value(&
            LDT_rc%lnc(n),LDT_rc%lnr(n),&
            CLSMF25_struc(n)%albvisdif%vlevels))
       
        ! Fill in derived parameter entries:
        ! ( input_parmattribs -> output_parmattribs ) 
       call populate_param_attribs( "ARS2CLSM", &
            "Catchment F2.5 wetness parameters","m2kg-1", &
            CLSMF25_struc(n)%ars1, &
            CLSMF25_struc(n)%ars2 )
       
       call populate_param_attribs( "ARS3CLSM", &
            "Catchment F2.5 wetness parameters","m2kg-1", &
            CLSMF25_struc(n)%ars1, &
            CLSMF25_struc(n)%ars3 )
       
       allocate(CLSMF25_struc(n)%ars1%value(&
            LDT_rc%lnc(n),LDT_rc%lnr(n),&
            CLSMF25_struc(n)%ars1%vlevels))
       allocate(CLSMF25_struc(n)%ars2%value(&
            LDT_rc%lnc(n),LDT_rc%lnr(n),&
            CLSMF25_struc(n)%ars2%vlevels))
       allocate(CLSMF25_struc(n)%ars3%value(&
            LDT_rc%lnc(n),LDT_rc%lnr(n),&
            CLSMF25_struc(n)%ars3%vlevels))
       
       ! Fill in derived parameter entries:
       ! ( input_parmattribs -> output_parmattribs ) 
       call populate_param_attribs( "ARA2CLSM", &
            "Catchment F2.5 topographic shape parameters","m2kg-1", &
            CLSMF25_struc(n)%ara1, &
            CLSMF25_struc(n)%ara2 )
       
       call populate_param_attribs( "ARA3CLSM", &
            "Catchment F2.5 topographic shape parameters","m2kg-1", &
            CLSMF25_struc(n)%ara1, &
            CLSMF25_struc(n)%ara3 )
       
       call populate_param_attribs( "ARA4CLSM", &
            "Catchment F2.5 topographic shape parameters","m2kg-1", &
            CLSMF25_struc(n)%ara1, &
            CLSMF25_struc(n)%ara4 )
       
       allocate(CLSMF25_struc(n)%ara1%value(&
            LDT_rc%lnc(n),LDT_rc%lnr(n),&
            CLSMF25_struc(n)%ara1%vlevels))
       allocate(CLSMF25_struc(n)%ara2%value(&
            LDT_rc%lnc(n),LDT_rc%lnr(n),&
            CLSMF25_struc(n)%ara2%vlevels))
       allocate(CLSMF25_struc(n)%ara3%value(&
            LDT_rc%lnc(n),LDT_rc%lnr(n),&
            CLSMF25_struc(n)%ara3%vlevels))
       allocate(CLSMF25_struc(n)%ara4%value(&
            LDT_rc%lnc(n),LDT_rc%lnr(n),&
            CLSMF25_struc(n)%ara4%vlevels))

       ! Fill in derived parameter entries:
       ! ( input_parmattribs -> output_parmattribs ) 
       call populate_param_attribs( "ARW2CLSM", &
            "Catchment F2.5 minimum theta parameters","m2kg-1", &
            CLSMF25_struc(n)%arw1, &
            CLSMF25_struc(n)%arw2 )
       
       call populate_param_attribs( "ARW3CLSM", &
            "Catchment F2.5 minimum theta parameters","m2kg-1", &
            CLSMF25_struc(n)%arw1, &
            CLSMF25_struc(n)%arw3 )
       
       call populate_param_attribs( "ARW4CLSM", &
            "Catchment F2.5 minimum theta parameters","m2kg-1", &
            CLSMF25_struc(n)%arw1, &
            CLSMF25_struc(n)%arw4 )          
       
       allocate(CLSMF25_struc(n)%arw1%value(&
            LDT_rc%lnc(n),LDT_rc%lnr(n),&
            CLSMF25_struc(n)%arw1%vlevels))
       allocate(CLSMF25_struc(n)%arw2%value(&
            LDT_rc%lnc(n),LDT_rc%lnr(n),&
            CLSMF25_struc(n)%arw2%vlevels))
       allocate(CLSMF25_struc(n)%arw3%value(&
            LDT_rc%lnc(n),LDT_rc%lnr(n),&
            CLSMF25_struc(n)%arw3%vlevels))
       allocate(CLSMF25_struc(n)%arw4%value(&
            LDT_rc%lnc(n),LDT_rc%lnr(n),&
            CLSMF25_struc(n)%arw4%vlevels))
       
       ! Fill in derived parameter entries:
       ! ( input_parmattribs -> output_parmattribs ) 
       call populate_param_attribs( "BF2CLSM", &
            "Catchment F2.5 baseflow topographic parameters","m", &
            CLSMF25_struc(n)%bf1, &
            CLSMF25_struc(n)%bf2 )
       
       call populate_param_attribs( "BF3CLSM", &
            "Catchment F2.5 baseflow topographic parameters","log(m)", &
            CLSMF25_struc(n)%bf1, &
            CLSMF25_struc(n)%bf3 )
       
       allocate(CLSMF25_struc(n)%bf1%value(&
            LDT_rc%lnc(n),LDT_rc%lnr(n),&
            CLSMF25_struc(n)%bf1%vlevels))
       allocate(CLSMF25_struc(n)%bf2%value(&
            LDT_rc%lnc(n),LDT_rc%lnr(n),&
            CLSMF25_struc(n)%bf2%vlevels))
       allocate(CLSMF25_struc(n)%bf3%value(&
            LDT_rc%lnc(n),LDT_rc%lnr(n),&
            CLSMF25_struc(n)%bf3%vlevels))
       ! Fill in derived parameter entries:
       ! ( input_parmattribs -> output_parmattribs ) 
       call populate_param_attribs( "TSA2CLSM", &
            "Catchment F2.5 water transfer parameters","-", &
            CLSMF25_struc(n)%tsa1, &
            CLSMF25_struc(n)%tsa2 )
       
       call populate_param_attribs( "TSB1CLSM", &
            "Catchment F2.5 water transfer parameters","-", &
            CLSMF25_struc(n)%tsa1, &
            CLSMF25_struc(n)%tsb1 )
       
       call populate_param_attribs( "TSB2CLSM", &
            "Catchment F2.5 water transfer parameters","-", &
            CLSMF25_struc(n)%tsa1, &
            CLSMF25_struc(n)%tsb2 )
       
       allocate(CLSMF25_struc(n)%tsa1%value(&
            LDT_rc%lnc(n),LDT_rc%lnr(n),&
            CLSMF25_struc(n)%tsa1%vlevels))
       allocate(CLSMF25_struc(n)%tsa2%value(&
            LDT_rc%lnc(n),LDT_rc%lnr(n),&
            CLSMF25_struc(n)%tsa2%vlevels))
       allocate(CLSMF25_struc(n)%tsb1%value(&
            LDT_rc%lnc(n),LDT_rc%lnr(n),&
            CLSMF25_struc(n)%tsb1%vlevels))
       allocate(CLSMF25_struc(n)%tsb2%value(&
            LDT_rc%lnc(n),LDT_rc%lnr(n),&
            CLSMF25_struc(n)%tsb2%vlevels))
       
       ! Fill in derived parameter entries:
       ! ( input_parmattribs -> output_parmattribs ) 
       call populate_param_attribs( "BTAUCLSM", &
            "Catchment F2.5 topographic tau parameters","-", &
            CLSMF25_struc(n)%atau, &
            CLSMF25_struc(n)%btau )
       
       allocate(CLSMF25_struc(n)%atau%value(&
            LDT_rc%lnc(n),LDT_rc%lnr(n),&
            CLSMF25_struc(n)%atau%vlevels))
       allocate(CLSMF25_struc(n)%btau%value(&
            LDT_rc%lnc(n),LDT_rc%lnr(n),&
            CLSMF25_struc(n)%btau%vlevels))
    enddo
 
!-- Catchment Parameter Config Entries:

!    call ESMF_ConfigFindLabel(LDT_config,"CLSMF25 tile coordinate file:",rc=rc)
!    do n=1,LDT_rc%nnest
!       call ESMF_ConfigGetAttribute(LDT_config,LDT_rc%tile_coord_file(n),rc=rc)
!       call LDT_verify(rc,'CLSMF25 tile coordinate file: not specified')
!    enddo

    call ESMF_ConfigFindLabel(LDT_config,"CLSMF25 soil param file:",rc=rc)
    do n=1,LDT_rc%nnest
       call ESMF_ConfigGetAttribute(LDT_config,&
            CLSMF25_struc(n)%soilparamfile,rc=rc)
       call LDT_verify(rc,'CLSMF25 soil param file: not specified')
    enddo

!    call ESMF_ConfigFindLabel(LDT_config,"CLSMF25 add to bedrock depth:",rc=rc)
!    do n=1,LDT_rc%nnest
!       call ESMF_ConfigGetAttribute(LDT_config,LDT_rc%addbdrckcrd(n),rc=rc)
!       call LDT_verify(rc,'CLSMF25 add to bedrock depth: not specified')
!       if( LDT_rc%addbdrckcrd(n) .ne. 0 ) then 
    CLSMF25_struc(:)%addbdrckcrd = 0.
!          write(LDT_logunit,*) " [INFO] ** The CLSM F2.5 Bedrock depth addition will "
!          write(LDT_logunit,*) " [INFO]  automatically be SET to 0. meters for now, "
!          write(LDT_logunit,*) " [INFO]  since this option may be removed at a later date. "
!       endif
!    enddo
    call ESMF_ConfigFindLabel(LDT_config,"CLSMF25 topo ar file:",rc=rc)
    do n=1,LDT_rc%nnest
       call ESMF_ConfigGetAttribute(LDT_config,CLSMF25_struc(n)%topo_ar_file,rc=rc)
       call LDT_verify(rc,'CLSMF25 topo ar file: not specified')
    enddo
    call ESMF_ConfigFindLabel(LDT_config,"CLSMF25 topo bf file:",rc=rc)
    do n=1,LDT_rc%nnest
       call ESMF_ConfigGetAttribute(LDT_config,CLSMF25_struc(n)%topo_bf_file,rc=rc)
       call LDT_verify(rc,'CLSMF25 topo bf file: not specified')
    enddo
    call ESMF_ConfigFindLabel(LDT_config,"CLSMF25 topo ts file:",rc=rc)
    do n=1,LDT_rc%nnest
       call ESMF_ConfigGetAttribute(LDT_config,CLSMF25_struc(n)%topo_ts_file,rc=rc)
       call LDT_verify(rc,'CLSMF25 topo ts file: not specified')
    enddo
    call ESMF_ConfigFindLabel(LDT_config,"CLSMF25 surf layer ts file:",rc=rc)
    do n=1,LDT_rc%nnest
       call ESMF_ConfigGetAttribute(LDT_config,CLSMF25_struc(n)%sltsfile,rc=rc)
       call LDT_verify(rc,'CLSMF25 surf layer ts file: not specified')
    enddo
    call ESMF_ConfigFindLabel(LDT_config,"CLSMF25 top soil layer depth:",rc=rc)
    do n=1,LDT_rc%nnest
       call ESMF_ConfigGetAttribute(LDT_config,CLSMF25_struc(n)%dzsfcrd,rc=rc)
       call LDT_verify(rc,'CLSMF25 top soil layer depth: not specified')
    enddo

    call ESMF_ConfigGetAttribute(LDT_config,catchparms_proj,&
         label="CLSMF25 map projection:",rc=rc)
    call LDT_verify(rc,'CLSMF25 map projection: option not specified in the config file')
    CLSMF25_struc(:)%catchparms_proj = catchparms_proj

    call ESMF_ConfigFindLabel(LDT_config,"CLSMF25 spatial transform:",rc=rc)
    do n=1,LDT_rc%nnest
       call ESMF_ConfigGetAttribute(LDT_config,&
            CLSMF25_struc(n)%catchparms_gridtransform,&
            rc=rc)
       call LDT_verify(rc,'CLSMF25 transform: option not specified in the config file')
    enddo
    
    call LDT_readDomainConfigSpecs("CLSMF25",catchparms_proj,&
         catchparms_gridDesc)
    
    do n=1,LDT_rc%nnest
       CLSMF25_struc(n)%catchparms_gridDesc = catchparms_gridDesc(n,:)
    enddo
!-- Read in Catchment Parameter Datasets:
    do n=1,LDT_rc%nnest

       call LDT_gridOptChecks( n, "CLSMF25", &
            CLSMF25_struc(n)%catchparms_gridtransform, &
            CLSMF25_struc(n)%catchparms_proj, &
            CLSMF25_struc(n)%catchparms_gridDesc(9) )
       
!       call LDT_catchparmsOptChecks( n, "CLSMF25", LDT_rc%catchparms_proj, &
!                    LDT_rc%catchparms_gridtransform )
       
    !- Tile coordinate data:
!       if( LDT_LSMparam_struc(n)%tilecoord%selectOpt == 1 ) then
!         write(LDT_logunit,*) 'Reading '//trim(LDT_rc%tile_coord_file(n))
!         call readtilecoord(trim(LDT_LSMparam_struc(n)%tilecoord%source)//char(0),&
!                          n,LDT_LSMparam_struc(n)%tilecoord%value(:,:,1))
!         write(LDT_logunit,*) 'Done reading '//trim(LDT_rc%tile_coord_file(n))
!       endif

    !- Vertical transmissivity (topo:ar) parameter data:
       write(LDT_logunit,*) "Reading gnu values from "//trim(CLSMF25_struc(n)%topo_ar_file)
       call read_CLSMF25_gnuparam(&
            n,CLSMF25_struc(n)%gnu%value(:,:,1),&
            LDT_LSMparam_struc(n)%landmask%value)
       write(LDT_logunit,*) "Done reading gnu values."
    !- Topographic parameters ("ar"):
       write(LDT_logunit,*) "Reading ars values from "//trim(CLSMF25_struc(n)%topo_ar_file)
       call read_CLSMF25_arsparams(&
            n,CLSMF25_struc(n)%ars1%value(:,:,1),&
            CLSMF25_struc(n)%ars2%value(:,:,1),&
            CLSMF25_struc(n)%ars3%value(:,:,1),&
            LDT_LSMparam_struc(n)%landmask%value)
       write(LDT_logunit,*) "Done reading ars values."
       write(LDT_logunit,*) "Reading ara values from "//trim(CLSMF25_struc(n)%topo_ar_file)
       call read_CLSMF25_araparams(&
            n,CLSMF25_struc(n)%ara1%value(:,:,1),&
            CLSMF25_struc(n)%ara2%value(:,:,1),&
            CLSMF25_struc(n)%ara3%value(:,:,1),&
            CLSMF25_struc(n)%ara4%value(:,:,1),&
            LDT_LSMparam_struc(n)%landmask%value)
       write(LDT_logunit,*) "Done reading ara values."
       write(LDT_logunit,*) "Reading arw values from "//trim(CLSMF25_struc(n)%topo_ar_file)
       call read_CLSMF25_arwparams(&
            n,CLSMF25_struc(n)%arw1%value(:,:,1),&
            CLSMF25_struc(n)%arw2%value(:,:,1),&
            CLSMF25_struc(n)%arw3%value(:,:,1),&
            CLSMF25_struc(n)%arw4%value(:,:,1),&
            LDT_LSMparam_struc(n)%landmask%value)
       write(LDT_logunit,*) "Done reading arw values."

    !- Baseflow timescale parameters:
       write(LDT_logunit,*) "Reading baseflow (bf) values from "//trim(CLSMF25_struc(n)%topo_bf_file)
       call read_CLSMF25_bfparams(&
            n,CLSMF25_struc(n)%bf1%value(:,:,1),&
            CLSMF25_struc(n)%bf2%value(:,:,1),&
            CLSMF25_struc(n)%bf3%value(:,:,1),&
            LDT_LSMparam_struc(n)%landmask%value)
       write(LDT_logunit,*) "Done reading baseflow (bf) values."

    !- Water transfer timescale parameters:
       write(LDT_logunit,*) "Reading ts values from "//trim(CLSMF25_struc(n)%topo_ts_file)
       call read_CLSMF25_tsparams(&
            n,CLSMF25_struc(n)%tsa1%value(:,:,1),&
            CLSMF25_struc(n)%tsa2%value(:,:,1),&
            CLSMF25_struc(n)%tsb1%value(:,:,1),&
            CLSMF25_struc(n)%tsb2%value(:,:,1),&
            LDT_LSMparam_struc(n)%landmask%value)
       write(LDT_logunit,*) "Done reading ts values."

    !- Surface layer timescale data:
       write(LDT_logunit,*) "Reading tau values from "//trim(CLSMF25_struc(n)%sltsfile)
       call read_CLSMF25_tauparams(&
            n,CLSMF25_struc(n)%atau%value(:,:,1),&
            CLSMF25_struc(n)%btau%value(:,:,1),&
            LDT_LSMparam_struc(n)%landmask%value)
       write(LDT_logunit,*) "Done reading tau values."


       call read_CLSMF25_psisat(&
            n,CLSMF25_struc(n)%psisat%value,&
            LDT_LSMparam_struc(n)%landmask%value)

       call read_CLSMF25_bexp(&
            n,CLSMF25_struc(n)%bexp%value,&
            LDT_LSMparam_struc(n)%landmask%value)

       call read_CLSMF25_wpwet(&
            n,CLSMF25_struc(n)%wpwet%value,&
            LDT_LSMparam_struc(n)%landmask%value)

       call read_CLSMF25_ksat(&
            n,CLSMF25_struc(n)%ksat%value,&
            LDT_LSMparam_struc(n)%landmask%value)
!       call read_CLSMF25_quartz(&
!            n,LDT_LSMparam_struc(n)%quartz%value,&
!            LDT_LSMparam_struc(n)%landmask%value)
!       call read_CLSMF25_soildepth(&
!            n,LDT_LSMparam_struc(n)%soildepth%value,&
!            LDT_LSMparam_struc(n)%landmask%value)
       
       call read_CLSMF25_bedrockdepth(&
            n,CLSMF25_struc(n)%bdrckdpth%value,&
            LDT_LSMparam_struc(n)%landmask%value)
       
       call populate_param_attribs( "ALBNIRDIFF", &
            "Alb near-IR diffuse scale factor", "-",  &
            CLSMF25_struc(n)%albnirdir, &
            CLSMF25_struc(n)%albnirdif )
       
       call populate_param_attribs( "ALBVISDIFF", &
            "Alb visible diffuse scale factor", "-",  &
            CLSMF25_struc(n)%albvisdir, &
            CLSMF25_struc(n)%albvisdif )

    enddo

    call ESMF_ConfigFindLabel(LDT_config,"Albedo NIR factor file:",rc=rc)
    do n=1,LDT_rc%nnest
       call ESMF_ConfigGetAttribute(LDT_config,&
            CLSMF25_struc(n)%albnirfile,rc=rc)
       call LDT_verify(rc,"Albedo NIR factor file: not specified")
    enddo

  ! Albedo VIS:
    call ESMF_ConfigFindLabel(LDT_config,"Albedo VIS factor file:",rc=rc)
    do n=1,LDT_rc%nnest
       call ESMF_ConfigGetAttribute(LDT_config,&
            CLSMF25_struc(n)%albvisfile,rc=rc)
       call LDT_verify(rc,"Albedo VIS factor file: not specified")
    enddo
    
    do n=1,LDT_rc%nnest
       write(LDT_logunit,*) 'Reading '//trim(CLSMF25_struc(n)%albnirfile)
       call read_CLSMF25_albnir(&
            n,CLSMF25_struc(n)%albnirdir%value, &   ! Direct
            CLSMF25_struc(n)%albnirdif%value, &   ! Diffuse
            LDT_LSMparam_struc(n)%landmask%value )
       write(LDT_logunit,*) 'Done reading '//&
            trim(CLSMF25_struc(n)%albnirfile)
       
       write(LDT_logunit,*) 'Reading '//trim(CLSMF25_struc(n)%albvisfile)
       call read_CLSMF25_albvis(&
            n,CLSMF25_struc(n)%albvisdir%value, &   ! Direct
            CLSMF25_struc(n)%albvisdif%value, &   ! Diffuse
            LDT_LSMparam_struc(n)%landmask%value )
       write(LDT_logunit,*) 'Done reading '//trim(CLSMF25_struc(n)%albvisfile)

       LDT_rc%monthlyData(n) = .true.
       LDT_albedo_struc(n)%albInterval = "monthly"
       LDT_gfrac_struc(n)%gfracInterval = "monthly"
    enddo

  end subroutine catchmentParms_init

  subroutine catchmentParms_writeHeader(n,ftn,dimID,monthID)

    integer     :: n
    integer     :: ftn
    integer     :: dimID(3)
    integer     :: monthID

    integer     :: t_dimID(3)
    integer     :: tdimID(3)

    tdimID(1) = dimID(1)
    tdimID(2) = dimID(2)

    t_dimID(1) = dimID(1)
    t_dimID(2) = dimID(2)

    if(LDT_gfrac_struc(n)%gfrac%selectOpt.gt.0) then
       if(LDT_gfrac_struc(n)%gfracInterval.eq."monthly") then !monthly
          t_dimID(3) = monthID
       endif
    end if

    call LDT_writeNETCDFdataHeader(n,ftn,dimID,&
         CLSMF25_struc(n)%gnu)
    call LDT_writeNETCDFdataHeader(n,ftn,dimID,&
         CLSMF25_struc(n)%ars1)
    call LDT_writeNETCDFdataHeader(n,ftn,dimID,&
         CLSMF25_struc(n)%ars2)
    call LDT_writeNETCDFdataHeader(n,ftn,dimID,&
         CLSMF25_struc(n)%ars3)
    call LDT_writeNETCDFdataHeader(n,ftn,dimID,&
         CLSMF25_struc(n)%ara1)
    call LDT_writeNETCDFdataHeader(n,ftn,dimID,&
         CLSMF25_struc(n)%ara2)
    call LDT_writeNETCDFdataHeader(n,ftn,dimID,&
         CLSMF25_struc(n)%ara3)
    call LDT_writeNETCDFdataHeader(n,ftn,dimID,&
         CLSMF25_struc(n)%ara4)
    
    call LDT_writeNETCDFdataHeader(n,ftn,dimID,&
         CLSMF25_struc(n)%arw1)
    call LDT_writeNETCDFdataHeader(n,ftn,dimID,&
         CLSMF25_struc(n)%arw2)
    call LDT_writeNETCDFdataHeader(n,ftn,dimID,&
         CLSMF25_struc(n)%arw3)
    call LDT_writeNETCDFdataHeader(n,ftn,dimID,&
         CLSMF25_struc(n)%arw4)
    
    call LDT_writeNETCDFdataHeader(n,ftn,dimID,&
         CLSMF25_struc(n)%bf1)
    call LDT_writeNETCDFdataHeader(n,ftn,dimID,&
         CLSMF25_struc(n)%bf2)
    call LDT_writeNETCDFdataHeader(n,ftn,dimID,&
         CLSMF25_struc(n)%bf3)
    
    call LDT_writeNETCDFdataHeader(n,ftn,dimID,&
         CLSMF25_struc(n)%tsa1)
    call LDT_writeNETCDFdataHeader(n,ftn,dimID,&
         CLSMF25_struc(n)%tsa2)
    call LDT_writeNETCDFdataHeader(n,ftn,dimID,&
         CLSMF25_struc(n)%tsb1)
    call LDT_writeNETCDFdataHeader(n,ftn,dimID,&
         CLSMF25_struc(n)%tsb2)

    call LDT_writeNETCDFdataHeader(n,ftn,dimID,&
         CLSMF25_struc(n)%atau)
    call LDT_writeNETCDFdataHeader(n,ftn,dimID,&
         CLSMF25_struc(n)%btau)

    call LDT_writeNETCDFdataHeader(n,ftn,tdimID,&
         CLSMF25_struc(n)%psisat)
    call LDT_writeNETCDFdataHeader(n,ftn,tdimID,&
         CLSMF25_struc(n)%ksat)
    call LDT_writeNETCDFdataHeader(n,ftn,tdimID,&
         CLSMF25_struc(n)%bexp)
    call LDT_writeNETCDFdataHeader(n,ftn,tdimID,&
         CLSMF25_struc(n)%wpwet)
!    call LDT_writeNETCDFdataHeader(n,ftn,tdimID,&
!         LDT_LSMparam_struc(n)%quartz)
!    call LDT_writeNETCDFdataHeader(n,ftn,tdimID,&
!         LDT_LSMparam_struc(n)%soildepth)
    call LDT_writeNETCDFdataHeader(n,ftn,tdimID,&
         CLSMF25_struc(n)%bdrckdpth)

 !- Albedo NIR scale factors:
    if( LDT_albedo_struc(n)%albInterval.eq."monthly" ) t_dimID(3) = monthID
    call LDT_writeNETCDFdataHeader(n,ftn,t_dimID,&
         CLSMF25_struc(n)%albnirdir)
    call LDT_writeNETCDFdataHeader(n,ftn,t_dimID,&
         CLSMF25_struc(n)%albnirdif)
    
 !- Albedo VIS scale factors:
    if( LDT_albedo_struc(n)%albInterval.eq."monthly" ) t_dimID(3) = monthID
    call LDT_writeNETCDFdataHeader(n,ftn,t_dimID,&
         CLSMF25_struc(n)%albvisdir)
    call LDT_writeNETCDFdataHeader(n,ftn,t_dimID,&
         CLSMF25_struc(n)%albvisdif)

    call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"ALBEDO_DATA_INTERVAL", &
         LDT_albedo_struc(n)%albInterval))
    
  end subroutine catchmentParms_writeHeader

  subroutine catchmentParms_writeData(n,ftn)

    integer   :: n 
    integer   :: ftn

    call LDT_writeNETCDFdata(n,ftn,CLSMF25_struc(n)%gnu)

    call LDT_writeNETCDFdata(n,ftn,CLSMF25_struc(n)%ars1)
    call LDT_writeNETCDFdata(n,ftn,CLSMF25_struc(n)%ars2)
    call LDT_writeNETCDFdata(n,ftn,CLSMF25_struc(n)%ars3)

    call LDT_writeNETCDFdata(n,ftn,CLSMF25_struc(n)%ara1)
    call LDT_writeNETCDFdata(n,ftn,CLSMF25_struc(n)%ara2)
    call LDT_writeNETCDFdata(n,ftn,CLSMF25_struc(n)%ara3)
    call LDT_writeNETCDFdata(n,ftn,CLSMF25_struc(n)%ara4)

    call LDT_writeNETCDFdata(n,ftn,CLSMF25_struc(n)%arw1)
    call LDT_writeNETCDFdata(n,ftn,CLSMF25_struc(n)%arw2)
    call LDT_writeNETCDFdata(n,ftn,CLSMF25_struc(n)%arw3)
    call LDT_writeNETCDFdata(n,ftn,CLSMF25_struc(n)%arw4)

    call LDT_writeNETCDFdata(n,ftn,CLSMF25_struc(n)%bf1)
    call LDT_writeNETCDFdata(n,ftn,CLSMF25_struc(n)%bf2)
    call LDT_writeNETCDFdata(n,ftn,CLSMF25_struc(n)%bf3)

    call LDT_writeNETCDFdata(n,ftn,CLSMF25_struc(n)%tsa1)
    call LDT_writeNETCDFdata(n,ftn,CLSMF25_struc(n)%tsa2)
    call LDT_writeNETCDFdata(n,ftn,CLSMF25_struc(n)%tsb1)
    call LDT_writeNETCDFdata(n,ftn,CLSMF25_struc(n)%tsb2)

    call LDT_writeNETCDFdata(n,ftn,CLSMF25_struc(n)%atau)
    call LDT_writeNETCDFdata(n,ftn,CLSMF25_struc(n)%btau)

    call LDT_writeNETCDFdata(n,ftn,CLSMF25_struc(n)%psisat)
    call LDT_writeNETCDFdata(n,ftn,CLSMF25_struc(n)%ksat)
    call LDT_writeNETCDFdata(n,ftn,CLSMF25_struc(n)%bexp)
    call LDT_writeNETCDFdata(n,ftn,CLSMF25_struc(n)%wpwet)
!    call LDT_writeNETCDFdata(n,ftn,LDT_LSMparam_struc(n)%quartz)
!    call LDT_writeNETCDFdata(n,ftn,LDT_LSMparam_struc(n)%soildepth)
    call LDT_writeNETCDFdata(n,ftn,CLSMF25_struc(n)%bdrckdpth)

    
    call LDT_writeNETCDFdata(n,ftn,CLSMF25_struc(n)%albnirdir)
    call LDT_writeNETCDFdata(n,ftn,CLSMF25_struc(n)%albnirdif)
    
    call LDT_writeNETCDFdata(n,ftn,CLSMF25_struc(n)%albvisdir)
    call LDT_writeNETCDFdata(n,ftn,CLSMF25_struc(n)%albvisdif)

  end subroutine catchmentParms_writeData


!BOP
! !ROUTINE:  set_param_attribs
! \label{set_param_attribs}
!
! !INTERFACE:
  subroutine set_param_attribs(paramEntry, short_name, vlevels, &
                 units, full_name )

! !DESCRIPTION:
!   This routine reads over the parameter attribute entries
!   in the param_attribs.txt file.
!
! !USES:
   type(LDT_paramEntry),intent(inout) :: paramEntry
   character(len=*),    intent(in)    :: short_name
   integer, optional                  :: vlevels
   character(len=*),     optional     :: units 
   character(len=*),     optional     :: full_name

   integer   :: v_temp
   character(20) :: unit_temp
   character(100):: name_temp

! ____________________________________________________
    
   if(present(vlevels)) then 
      v_temp = vlevels
   else
      v_temp = 1
   endif

   if(present(units)) then
      unit_temp = units
   else
      unit_temp = "none"
   endif

   if(present(full_name)) then
      name_temp = full_name
   else
      name_temp = trim(short_name)
   endif

   paramEntry%short_name = trim(short_name)
   paramEntry%vlevels = v_temp
   paramEntry%selectOpt = 1
   paramEntry%source = "CLSMF2.5"
   paramEntry%units = trim(unit_temp)
   paramEntry%num_times = 1
   paramEntry%num_bins = 1
   paramEntry%standard_name = trim(name_temp)

  end subroutine set_param_attribs

end module CLSMF25_parmsMod
