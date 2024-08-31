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
module JULES50_parmsMod
!BOP
!
! !MODULE: JULES50_parmsMod
!
! !DESCRIPTION:
!  The code in this file implements routines to read JULES50 parameter
!  data. 
!  \subsubsection{Overview}
!  This routines in this module provides routines to read the 
!  JULES50 parameter file data.
!
! !REVISION HISTORY:
!
!  04 Nov 2013: K. Arsenault: Added layers for JULES50 model
!  06 Apr 2017: Shugong Wang: Add soil parameters for JULES model 
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif

  use ESMF
  use LDT_coreMod
  use LDT_historyMod
  use LDT_paramDataMod
  use ldt_logmod
  use LDT_paramMaskCheckMod

  implicit none

  PRIVATE
!------------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!------------------------------------------------------------------------------
  public :: JULES50parms_init    
  public :: JULES50parms_writeHeader
  public :: JULES50parms_writeData

!------------------------------------------------------------------------------
! !PUBLIC TYPES:
!------------------------------------------------------------------------------
  public :: JULES50_struc

  type, public :: jules_type_dec
     real      :: jules_undef

     real           :: soiltype_gridDesc(20)
     character*50   :: soiltype_gridtransform
     character*50   :: soilalbedo_gridtransform
! -  JULES50-specific:
     type(LDT_paramEntry) :: jules  ! JULES50 parameters (collective)
     type(LDT_paramEntry) :: b_map 
     type(LDT_paramEntry) :: sathh_map 
     type(LDT_paramEntry) :: satcon_map 
     type(LDT_paramEntry) :: sm_sat_map 
     type(LDT_paramEntry) :: sm_wilt_map 
     type(LDT_paramEntry) :: sm_crit_map 
     type(LDT_paramEntry) :: hcap_map 
     type(LDT_paramEntry) :: hcon_map
     type(LDT_paramEntry) :: usda_type
     type(LDT_paramEntry) :: soil_albedo
     type(LDT_paramEntry) :: top_ti_mean
     type(LDT_paramEntry) :: top_ti_sig
     type(LDT_paramEntry) :: top_fexp

     real, allocatable :: soil_type(:,:) 
     character(len=128)   :: nc_cap_param    ! NetCDF file of CAP soil texture fraction file 
     character(len=128)   :: nc_whs_param    ! NetCDF file of WHS soil albedo input file 
     !character(len=32)    :: soil_param_mode
     character(len=128)   :: nc_um_ancillary
  end type jules_type_dec

  type(jules_type_dec), allocatable :: JULES50_struc(:)

  type(LDT_fillopts) :: usda_fill_opts
  type(LDT_fillopts) :: albedo_fill_opts
!  real, private :: organic_soil(8,3)  = reshape((/2.7, 6.1, 12.0,             & ! b
!                                                  0.0103, 0.0102, 0.0101,     & ! sathh
!                                                  0.28, 0.002, 0.0001,        & ! satcon 
!                                                  0.93, 0.88, 0.83,           & ! sm_sat
!                                                  0.11, 0.34, 0.51,           & ! sm_crit
!                                                  0.03, 0.18, 0.37,           & ! sm_wilt
!                                                  0.58e+6, 0.58e+6, 0.58e+6,  & ! hcap
!                                                  0.06, 0.06, 0.06/),         & ! hcon
!                                                  (/8, 3/))
contains

!BOP
! 
! !ROUTINE: JULES50parms_init
! \label{JULES50parms_init}
! 
! !INTERFACE:
  subroutine JULES50parms_init(flag)

! !USES:
   use LDT_logMod,    only : LDT_verify, LDT_endrun, &
             LDT_getNextUnitNumber, LDT_releaseUnitNumber
  use LDT_coreMod
!
! !DESCRIPTION:
!
! Joint UK Land Environment Simulator, v5.0 (JULES50) model parameters.
! 
!  The routines invoked are: 
!  \begin{description}
!   \item[julesParmssetup](\ref{julesParmssetup}) \newline
!    calls the registry to invoke the julesParms setup methods. 
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
   character(len=32)    :: soil_param_mode

   ! _____________________________________________________________________

   allocate(JULES50_struc(LDT_rc%nnest))

   do n=1,LDT_rc%nnest
      ! - JULES50 parameters:
      call set_param_attribs(JULES50_struc(n)%jules,"JULES50")
      call set_param_attribs(JULES50_struc(n)%b_map      ,"b")
      call set_param_attribs(JULES50_struc(n)%sathh_map  ,"sathh")
      call set_param_attribs(JULES50_struc(n)%satcon_map ,"satcon")
      call set_param_attribs(JULES50_struc(n)%sm_sat_map ,"sm_sat")
      call set_param_attribs(JULES50_struc(n)%sm_wilt_map,"sm_wilt")
      call set_param_attribs(JULES50_struc(n)%sm_crit_map,"sm_crit")
      call set_param_attribs(JULES50_struc(n)%hcap_map   ,"hcap")
      call set_param_attribs(JULES50_struc(n)%hcon_map   ,"hcon")
      call set_param_attribs(JULES50_struc(n)%usda_type  ,"usda_type")
      call set_param_attribs(JULES50_struc(n)%soil_albedo,"soil_albedo")      
      call set_param_attribs(JULES50_struc(n)%top_ti_mean,"top_ti_mean")
      call set_param_attribs(JULES50_struc(n)%top_ti_sig ,"top_ti_sig")
      call set_param_attribs(JULES50_struc(n)%top_fexp   ,"top_fexp")

      allocate(JULES50_struc(n)%soil_type(LDT_rc%lnc(n),LDT_rc%lnr(n)))
      call populate_param_attribs( "JULES_B",       "b parameter"," -", JULES50_struc(n)%jules, JULES50_struc(n)%b_map )
      call populate_param_attribs( "JULES_SATHH",   "sathh",      " -", JULES50_struc(n)%jules, JULES50_struc(n)%sathh_map )
      call populate_param_attribs( "JULES_SATCON",  "satcon",     " -", JULES50_struc(n)%jules, JULES50_struc(n)%satcon_map )
      call populate_param_attribs( "JULES_SM_SAT",  "sm_sat",     " -", JULES50_struc(n)%jules, JULES50_struc(n)%sm_sat_map )
      call populate_param_attribs( "JULES_SM_CRIT", "sm_crit",    " -", JULES50_struc(n)%jules, JULES50_struc(n)%sm_crit_map )
      call populate_param_attribs( "JULES_SM_WILT", "sm_wilt",    " -", JULES50_struc(n)%jules, JULES50_struc(n)%sm_wilt_map )
      call populate_param_attribs( "JULES_HCAP",    "hcap",       " -", JULES50_struc(n)%jules, JULES50_struc(n)%hcap_map )
      call populate_param_attribs( "JULES_HCON",    "hcon",       " -", JULES50_struc(n)%jules, JULES50_struc(n)%hcon_map )
      call populate_param_attribs( "USDA_SOIL_TYPE","soil type",  " -", JULES50_struc(n)%jules, JULES50_struc(n)%usda_type )
      call populate_param_attribs( "JULES_ALBSOIL",  "soil albedo"," -", JULES50_struc(n)%jules, JULES50_struc(n)%soil_albedo )
      call populate_param_attribs( "JULES_ALBSOIL", "soil albedo"," -", JULES50_struc(n)%jules, JULES50_struc(n)%soil_albedo )
      call populate_param_attribs( "JULES_TI_MEAN", "TI mean",    " -", JULES50_struc(n)%jules, JULES50_struc(n)%top_ti_mean )
      call populate_param_attribs( "JULES_TI_SIG",  "TI std",     " -", JULES50_struc(n)%jules, JULES50_struc(n)%top_ti_sig )
      call populate_param_attribs( "JULES_FEXP",    "fexp",       " -", JULES50_struc(n)%jules, JULES50_struc(n)%top_fexp )

      allocate(JULES50_struc(n)%b_map%value(LDT_rc%lnc(n),LDT_rc%lnr(n),JULES50_struc(n)%b_map%vlevels))       
      allocate(JULES50_struc(n)%sathh_map%value(LDT_rc%lnc(n),LDT_rc%lnr(n),JULES50_struc(n)%sathh_map%vlevels))       
      allocate(JULES50_struc(n)%satcon_map%value(LDT_rc%lnc(n),LDT_rc%lnr(n),JULES50_struc(n)%satcon_map%vlevels))       
      allocate(JULES50_struc(n)%sm_sat_map%value(LDT_rc%lnc(n),LDT_rc%lnr(n),JULES50_struc(n)%sm_sat_map%vlevels))       
      allocate(JULES50_struc(n)%sm_wilt_map%value(LDT_rc%lnc(n),LDT_rc%lnr(n),JULES50_struc(n)%sm_wilt_map%vlevels))       
      allocate(JULES50_struc(n)%sm_crit_map%value(LDT_rc%lnc(n),LDT_rc%lnr(n),JULES50_struc(n)%sm_crit_map%vlevels))       
      allocate(JULES50_struc(n)%hcap_map%value(LDT_rc%lnc(n),LDT_rc%lnr(n),JULES50_struc(n)%hcap_map%vlevels))       
      allocate(JULES50_struc(n)%hcon_map%value(LDT_rc%lnc(n),LDT_rc%lnr(n),JULES50_struc(n)%hcon_map%vlevels))       
      allocate(JULES50_struc(n)%usda_type%value(LDT_rc%lnc(n),LDT_rc%lnr(n),JULES50_struc(n)%usda_type%vlevels))       
      allocate(JULES50_struc(n)%soil_albedo%value(LDT_rc%lnc(n),LDT_rc%lnr(n),JULES50_struc(n)%soil_albedo%vlevels))       
      allocate(JULES50_struc(n)%top_ti_mean%value(LDT_rc%lnc(n),LDT_rc%lnr(n),JULES50_struc(n)%top_ti_mean%vlevels))       
      allocate(JULES50_struc(n)%top_ti_sig%value(LDT_rc%lnc(n),LDT_rc%lnr(n),JULES50_struc(n)%top_ti_sig%vlevels))       
      allocate(JULES50_struc(n)%top_fexp%value(LDT_rc%lnc(n),LDT_rc%lnr(n),JULES50_struc(n)%top_fexp%vlevels))   
   enddo

   write(LDT_logunit,*)" - - - - - - - - - - JULES v5.0 LSM Parameters - - - - - - - - - - - - -"
    
   call ESMF_ConfigFindLabel(LDT_config,"JULES soil parameter mode:",rc=rc)
   do n=1,LDT_rc%nnest
      call ESMF_ConfigGetAttribute(LDT_config, soil_param_mode,rc=rc)
      call LDT_verify(rc,'JULES soil parameter mode: not specified')
   enddo

   if(trim(soil_param_mode) .eq. "readin") then
     call ESMF_ConfigFindLabel(LDT_config,"JULES ancillary file:",rc=rc)
     do n=1,LDT_rc%nnest
        call ESMF_ConfigGetAttribute(LDT_config, JULES50_struc(n)%nc_um_ancillary,rc=rc)
        call LDT_verify(rc,'JULES ancillary file: not specified')

        inquire(file=trim(JULES50_struc(n)%nc_um_ancillary), exist=file_exists)
        if( .not. file_exists ) then
           write(LDT_logunit,*) "[ERR] JULES ancillary file ",&
                trim(JULES50_struc(n)%nc_um_ancillary)," does not exist."
           call LDT_endrun
        endif
     enddo
   else
     call ESMF_ConfigFindLabel(LDT_config,"JULES WHS soil parameter file:",rc=rc)
     do n=1,LDT_rc%nnest
        call ESMF_ConfigGetAttribute(LDT_config, JULES50_struc(n)%nc_whs_param,rc=rc)
        call LDT_verify(rc,'JULES WHS soil parameter file: not specified')

        inquire(file=trim(JULES50_struc(n)%nc_whs_param), exist=file_exists)
        if( .not. file_exists ) then
           write(LDT_logunit,*) "[ERR] JULES WHS soil parameter file ",&
                trim(JULES50_struc(n)%nc_whs_param)," does not exist."
           call LDT_endrun
        endif
     enddo
     
     call ESMF_ConfigFindLabel(LDT_config,"JULES CAP soil fraction file:",rc=rc)
     do n=1,LDT_rc%nnest
        call ESMF_ConfigGetAttribute(LDT_config, JULES50_struc(n)%nc_cap_param,rc=rc)
        call LDT_verify(rc,'JULES CAP soil fraction file: not specified')

        inquire(file=trim(JULES50_struc(n)%nc_cap_param), exist=file_exists)
        if( .not. file_exists ) then
           write(LDT_logunit,*) "[ERR] JULES CAP Soil Fraction File ",&
                trim(JULES50_struc(n)%nc_cap_param)," does not exist."
           call LDT_endrun
        endif
     enddo
     
     ! Read in soil type "fill" options:
     usda_fill_opts%filltype = "none"
     call ESMF_ConfigGetAttribute(LDT_config, usda_fill_opts%filltype, & 
          label="Soil texture fill option:",rc=rc)
     call LDT_verify(rc,"Soil texture fill option: option not specified in the config file")
     
     if( usda_fill_opts%filltype == "neighbor" ) then
        call ESMF_ConfigGetAttribute(LDT_config, usda_fill_opts%fillradius, &
             label="Soil texture fill radius:",rc=rc)
        call LDT_verify(rc,"Soil texture fill radius: option not specified in the config file")
        
        call ESMF_ConfigGetAttribute(LDT_config, usda_fill_opts%fillvalue, &
             label="Soil texture fill value:",rc=rc)
        call LDT_verify(rc,"Soil texture fill value: option not specified in the config file")
        if( usda_fill_opts%fillvalue > 12. ) then
           usda_fill_opts%fillvalue = 5.  !silty loam  
        end if
     elseif( usda_fill_opts%filltype == "none" ) then
        write(LDT_logunit,*) "[INFO] 'NONE' Parameter-Mask Agreement Option Selected for Soil texture"
     else
        write(LDT_logunit,*) "[ERR] Fill option for soil texture is not valid: ",trim(usda_fill_opts%filltype)
        write(LDT_logunit,*) "  Please select one of these:  none or neighbor "
        write(LDT_logunit,*) "  Programming stopping ..."
        call LDT_endrun
     end if
    
   ! Read in soil albedo "fill" options:
     albedo_fill_opts%filltype = "none"
     call ESMF_ConfigGetAttribute(LDT_config, albedo_fill_opts%filltype, & 
          label="Albedo fill option:",rc=rc)
     call LDT_verify(rc,"Albedo fill option: option not specified in the config file")
     
     if( albedo_fill_opts%filltype == "neighbor" ) then
        call ESMF_ConfigGetAttribute(LDT_config, albedo_fill_opts%fillradius, &
             label="Albedo fill radius:",rc=rc)
        call LDT_verify(rc,"Albedo fill radius: option not specified in the config file")
        
        call ESMF_ConfigGetAttribute(LDT_config, albedo_fill_opts%fillvalue, &
             label="Albedo fill value:",rc=rc)
        call LDT_verify(rc,"Albedo fill value: option not specified in the config file")
        if( albedo_fill_opts%fillvalue > 1.0) then
           albedo_fill_opts%fillvalue = 0.14  
        end if
     elseif( albedo_fill_opts%filltype == "none" ) then
        write(LDT_logunit,*) "[INFO] 'NONE' Parameter-Mask Agreement Option Selected for Albedo"
     else
        write(LDT_logunit,*) "[ERR] Fill option for albedo is not valid: ",trim(albedo_fill_opts%filltype)
        write(LDT_logunit,*) "  Please select one of these:  none or neighbor "
        write(LDT_logunit,*) "  Programming stopping ..."
        call LDT_endrun
     end if

     call ESMF_ConfigFindLabel(LDT_config,"Soil texture spatial transform:",rc=rc)
     do n=1,LDT_rc%nnest
        call ESMF_ConfigGetAttribute(LDT_config, JULES50_struc(n)%soiltype_gridtransform, &
             rc=rc)
        call LDT_verify(rc,'Soil texture spatial transform: option not specified in the config file')
     enddo
     
     call ESMF_ConfigFindLabel(LDT_config,"Albedo spatial transform:",rc=rc)
     do n=1,LDT_rc%nnest
        call ESMF_ConfigGetAttribute(LDT_config, JULES50_struc(n)%soilalbedo_gridtransform, &
             rc=rc)
        call LDT_verify(rc,'Albedo spatial transform: option not specified in the config file')
     enddo
     write(LDT_logunit,*)"[INFO] ..."
   endif


   if(trim(soil_param_mode) .eq. "readin") then
     call JULES50_read_soil_params 
   else 
     call JULES50_calc_soil_params
   endif

 end subroutine JULES50parms_init
 
 subroutine JULES50_read_soil_params
  use LDT_logMod
  use LDT_paramDataMod
  use LDT_coreMod
  
  integer :: n, col, row, soil_id
  real, allocatable :: param_array(:,:)
   
  do n=1,LDT_rc%nnest
    allocate(param_array(LDT_rc%lnc(n),LDT_rc%lnr(n)))
    ! firstly assign soil_type to -9999
    do col=1, LDT_rc%lnc(n)
      do row=1, LDT_rc%lnr(n)
        JULES50_struc(n)%usda_type%value(col, row, 1)   = -9999 
      enddo
    enddo
    
    call read_um_soil_params(n, trim(JULES50_struc(n)%nc_um_ancillary), "albedo", param_array) 
    JULES50_struc(n)%soil_albedo%value(:, :, 1)       = param_array(:,:) 
    
    call read_um_soil_params(n, trim(JULES50_struc(n)%nc_um_ancillary), "smvcst", param_array) 
    JULES50_struc(n)%sm_sat_map%value(:, :, 1)       = param_array(:,:) 
    
    call read_um_soil_params(n, trim(JULES50_struc(n)%nc_um_ancillary), "smvccl", param_array) 
    JULES50_struc(n)%sm_crit_map%value(:, :, 1)       = param_array(:,:) 
    
    call read_um_soil_params(n, trim(JULES50_struc(n)%nc_um_ancillary), "smvcwt", param_array) 
    JULES50_struc(n)%sm_wilt_map%value(:, :, 1)       = param_array(:,:) 
    
    call read_um_soil_params(n, trim(JULES50_struc(n)%nc_um_ancillary), "satcon", param_array) 
    JULES50_struc(n)%satcon_map%value(:, :, 1)       = param_array(:,:) 
    
    call read_um_soil_params(n, trim(JULES50_struc(n)%nc_um_ancillary), "sathh", param_array) 
    JULES50_struc(n)%sathh_map%value(:, :, 1)       = param_array(:,:) 
    
    call read_um_soil_params(n, trim(JULES50_struc(n)%nc_um_ancillary), "hcap", param_array) 
    JULES50_struc(n)%hcap_map%value(:, :, 1)       = param_array(:,:) 
    
    call read_um_soil_params(n, trim(JULES50_struc(n)%nc_um_ancillary), "hcon", param_array) 
    JULES50_struc(n)%hcon_map%value(:, :, 1)       = param_array(:,:) 
    
    ! b parameter 
    call read_um_soil_params(n, trim(JULES50_struc(n)%nc_um_ancillary), "b", param_array) 
    JULES50_struc(n)%b_map%value(:, :, 1)       = param_array(:,:) 

    ! topmodel
    call read_um_soil_params(n, trim(JULES50_struc(n)%nc_um_ancillary), "mean_ti", param_array) 
    JULES50_struc(n)%top_ti_mean%value(:, :, 1)       = param_array(:,:) 
  
    call read_um_soil_params(n, trim(JULES50_struc(n)%nc_um_ancillary), "std_ti", param_array) 
    JULES50_struc(n)%top_ti_sig%value(:, :, 1)       = param_array(:,:) 
    
    call read_um_soil_params(n, trim(JULES50_struc(n)%nc_um_ancillary), "fexp", param_array) 
    JULES50_struc(n)%top_fexp%value(:, :, 1)       = param_array(:,:) 
    
    deallocate(param_array) 
 enddo

 end subroutine JULES50_read_soil_params


 subroutine JULES50_calc_soil_params
  use LDT_logMod
  use LDT_paramDataMod
  use LDT_coreMod
  integer :: n, col, row, soil_id
  real :: b, sathh, satcon, sm_sat, sm_wilt, sm_crit, hcap, hcon
  do n=1,LDT_rc%nnest
    call create_soil_params_from_cap(n, JULES50_struc(n)%nc_cap_param, JULES50_struc(n)%soil_type)
    call create_soil_albedo_from_whs(n, JULES50_struc(n)%nc_whs_param, JULES50_struc(n)%soil_albedo%value)
    ! firstly assign soil_type to usda_type
    do col=1, LDT_rc%lnc(n)
      do row=1, LDT_rc%lnr(n)
        JULES50_struc(n)%usda_type%value(col, row, 1)   = JULES50_struc(n)%soil_type(col, row)
      enddo
    enddo
        
    !Fill where soil albedo values are missing compared to land/water mask:
    if( albedo_fill_opts%filltype == "neighbor" .or. &
        albedo_fill_opts%filltype == "average" ) then
      write(LDT_logunit,*) "Checking/filling mask values for: ", &
                        trim(JULES50_struc(n)%soil_albedo%short_name)
      write(fill_logunit,*) "Checking/filling mask values for: ", &
                        trim(JULES50_struc(n)%soil_albedo%short_name)
      albedo_fill_opts%watervalue = LDT_rc%udef
      call LDT_contIndivParam_Fill( n, LDT_rc%lnc(n), LDT_rc%lnr(n), &
            JULES50_struc(n)%soilalbedo_gridtransform,                   &
            JULES50_struc(n)%soil_albedo%vlevels,              &
            JULES50_struc(n)%soil_albedo%value, albedo_fill_opts%watervalue, &
            LDT_LSMparam_struc(n)%landmask2%value,              &
            albedo_fill_opts%filltype, albedo_fill_opts%fillvalue, albedo_fill_opts%fillradius )
    endif
    
    ! then, fill usda_type
    ! Fill where parameter values are missing compared to land/water mask:
    if( usda_fill_opts%filltype == "neighbor" ) then
       write(LDT_logunit,*) "Checking/filling mask values for: ", &
            trim(JULES50_struc(n)%usda_type%short_name)
       write(fill_logunit,*) "Checking/filling mask values for: ", &
            trim(JULES50_struc(n)%usda_type%short_name)
       !usda_fill_opts%watervalue = 7.
       usda_fill_opts%watervalue = -9999 ! this is the water type defined in soil category 
       call LDT_discreteParam_Fill( n, LDT_rc%lnc(n), LDT_rc%lnr(n), &
            JULES50_struc(n)%soiltype_gridtransform,                   &
            JULES50_struc(n)%usda_type%num_bins,                    &
            JULES50_struc(n)%usda_type%value, usda_fill_opts%watervalue, &
            LDT_LSMparam_struc(n)%landmask2%value,               &
            usda_fill_opts%filltype, usda_fill_opts%fillvalue, usda_fill_opts%fillradius )
    endif

    do col=1, LDT_rc%lnc(n)
      do row=1, LDT_rc%lnr(n)
        !soil_id = int(JULES50_struc(n)%soil_type(col, row))
        soil_id = int(JULES50_struc(n)%usda_type%value(col, row, 1))
        call calc_soil_params_id(soil_id, b, sathh, satcon, sm_sat, sm_wilt, sm_crit, hcap, hcon)
        JULES50_struc(n)%b_map%value(col, row, 1)       = b 
        JULES50_struc(n)%sathh_map%value(col, row, 1)   = sathh 
        JULES50_struc(n)%satcon_map%value(col, row, 1)  = satcon 
        JULES50_struc(n)%sm_sat_map%value(col, row, 1)  = sm_sat
        JULES50_struc(n)%sm_wilt_map%value(col, row, 1) = sm_wilt 
        JULES50_struc(n)%sm_crit_map%value(col, row, 1) = sm_crit
        JULES50_struc(n)%hcap_map%value(col, row, 1)    = hcap 
        JULES50_struc(n)%hcon_map%value(col, row, 1)    = hcon
        !JULES50_struc(n)%usda_type%value(col, row, 1)   = JULES50_struc(n)%soil_type(col, row)
      
        !!!! set soil parameters to 0 for land ice grid boxes
        if(LDT_LSMparam_struc(n)%landcover%value(col,row,9) .eq. 1.0) then
          JULES50_struc(n)%b_map%value(col, row, 1)       = 0.0
          JULES50_struc(n)%sathh_map%value(col, row, 1)   = 0.0
          JULES50_struc(n)%satcon_map%value(col, row, 1)  = 0.0
          JULES50_struc(n)%sm_sat_map%value(col, row, 1)  = 0.0
          JULES50_struc(n)%sm_wilt_map%value(col, row, 1) = 0.0
          JULES50_struc(n)%sm_crit_map%value(col, row, 1) = 0.0
          JULES50_struc(n)%hcap_map%value(col, row, 1)    = 630000.0 ! unit? 
          JULES50_struc(n)%hcon_map%value(col, row, 1)    = 0.265    ! uint? 
          JULES50_struc(n)%soil_albedo%value(col, row, 1) = 0.75     !
        endif
      enddo
    enddo
  enddo
 end subroutine JULES50_calc_soil_params

 subroutine JULES50parms_writeHeader(n, ftn, dimID, monthID)
   
   integer   :: n 
   integer   :: ftn
   integer   :: dimID(3)
   integer   :: monthID
   
   integer   :: t_dimID(3)
   
   if( JULES50_struc(n)%jules%selectOpt == 1 ) then
      call LDT_writeNETCDFdataHeader(n, ftn, dimID, JULES50_struc(n)%b_map)
      call LDT_writeNETCDFdataHeader(n, ftn, dimID, JULES50_struc(n)%sathh_map)
      call LDT_writeNETCDFdataHeader(n, ftn, dimID, JULES50_struc(n)%satcon_map)
      call LDT_writeNETCDFdataHeader(n, ftn, dimID, JULES50_struc(n)%sm_sat_map)
      call LDT_writeNETCDFdataHeader(n, ftn, dimID, JULES50_struc(n)%sm_wilt_map)
      call LDT_writeNETCDFdataHeader(n, ftn, dimID, JULES50_struc(n)%sm_crit_map)
      call LDT_writeNETCDFdataHeader(n, ftn, dimID, JULES50_struc(n)%hcap_map)
      call LDT_writeNETCDFdataHeader(n, ftn, dimID, JULES50_struc(n)%hcon_map)
      call LDT_writeNETCDFdataHeader(n, ftn, dimID, JULES50_struc(n)%usda_type)
      call LDT_writeNETCDFdataHeader(n, ftn, dimID, JULES50_struc(n)%soil_albedo)
      call LDT_writeNETCDFdataHeader(n, ftn, dimID, JULES50_struc(n)%top_ti_mean)
      call LDT_writeNETCDFdataHeader(n, ftn, dimID, JULES50_struc(n)%top_ti_sig)
      call LDT_writeNETCDFdataHeader(n, ftn, dimID, JULES50_struc(n)%top_fexp)
   endif
 end subroutine JULES50parms_writeHeader
 
 subroutine JULES50parms_writeData(n,ftn)

   integer   :: n 
   integer   :: ftn

   if( JULES50_struc(n)%jules%selectOpt == 1 ) then
      call LDT_writeNETCDFdata(n, ftn, JULES50_struc(n)%b_map)
      call LDT_writeNETCDFdata(n, ftn, JULES50_struc(n)%sathh_map)
      call LDT_writeNETCDFdata(n, ftn, JULES50_struc(n)%satcon_map)
      call LDT_writeNETCDFdata(n, ftn, JULES50_struc(n)%sm_sat_map)
      call LDT_writeNETCDFdata(n, ftn, JULES50_struc(n)%sm_wilt_map)
      call LDT_writeNETCDFdata(n, ftn, JULES50_struc(n)%sm_crit_map)
      call LDT_writeNETCDFdata(n, ftn, JULES50_struc(n)%hcap_map)
      call LDT_writeNETCDFdata(n, ftn, JULES50_struc(n)%hcon_map)
      call LDT_writeNETCDFdata(n, ftn, JULES50_struc(n)%usda_type)
      call LDT_writeNETCDFdata(n, ftn, JULES50_struc(n)%soil_albedo)
      call LDT_writeNETCDFdata(n, ftn, JULES50_struc(n)%top_ti_mean)
      call LDT_writeNETCDFdata(n, ftn, JULES50_struc(n)%top_ti_sig)
      call LDT_writeNETCDFdata(n, ftn, JULES50_struc(n)%top_fexp)

   endif

 end subroutine JULES50parms_writeData

!BOP
! !ROUTINE:  set_param_attribs
! \label{set_param_attribs}
!
! !INTERFACE:
  subroutine set_param_attribs(paramEntry, short_name)

! !DESCRIPTION:
!   This routine reads over the parameter attribute entries
!   in the param_attribs.txt file.
!
! !USES:
    type(LDT_paramEntry),intent(inout) :: paramEntry
    character(len=*),    intent(in)    :: short_name

! ____________________________________________________
    
   
    paramEntry%short_name = trim(short_name)
    paramEntry%vlevels = 1
    paramEntry%selectOpt = 1
    paramEntry%source = "JULES50"
    paramEntry%units ="none"
    paramEntry%num_times = 1
    paramEntry%num_bins = 1
    paramEntry%standard_name = trim(short_name)

  end subroutine set_param_attribs
    
  ! Calculate soil parameters for non-organic soil
  ! using [AncilDoc_Cosby], which references [Cosby 1984]
  subroutine calc_soil_params_id(soil_id, b, sathh, satcon, sm_sat, sm_wilt, sm_crit, hcap, hcon)
    implicit none
    integer :: soil_id 
    real :: b, sathh, satcon, sm_sat, sm_wilt, sm_crit, hcap, hcon
    real :: lambda_air, lambda_clay, lambda_sand, lambda_silt, aa, bb, cc
    real :: c_clay, c_sand, c_silt, f_silt, f_sand, f_clay 
  
    real :: soil_texture(12, 3) 
    
    soil_texture = transpose(reshape((/ 5, 92,  3, &       !  'Sand':           
                                       12, 82,  6, &       !  'Loamy sand':     
                                       32, 58, 10, &       !  'Sandy loam':     
                                       39, 43, 18, &       !  'Loam':           
                                       70, 17, 13, &       !  'Silty loam':     
                                       15, 58, 27, &       !  'Sandy clay loam':
                                       34, 32, 34, &       !  'Clay loam':      
                                       56, 10, 34, &       !  'Silty clay loam':
                                        6, 52, 42, &       !  'Sandy clay':     
                                       47,  6, 47, &       !  'Silty clay':     
                                       20, 22, 58, &       !  'Clay':           
                                       90,  7,  3 /), &    !  'Silt':           
                                       (/3, 12/)))
    
    if(soil_id>0) then
      f_silt = soil_texture(soil_id,1)*0.01;
      f_sand = soil_texture(soil_id,2)*0.01;
      f_clay = soil_texture(soil_id,3)*0.01;
      
      ! Clapp Hornberger parameter b
      ! JULES_SOIL_PROPS units: dimensionless.
      b      =   3.10 + 15.70 * f_clay - 0.3 * f_sand

      ! Saturated soil water suction SATHH in terms of log to the base 10
      ! (n.b. used to use nat log)
      ! JULES_SOIL_PROPS units: m. (says 'the *absolute* value of the soil
      !  matrix suction at saturation')
      sathh = 0.01 * (10.0 ** (2.17 - 0.63 * f_clay - 1.58 * f_sand))
        
      ! Saturated hydrological soil conductivity Ks in terms of log to the
      ! base 10 (n.b. used to use nat log)
      ! JULES_SOIL_PROPS units: kg m^-2 s^-1
      satcon = 10.0 ** (-2.75 - 0.64 * f_clay + 1.26 * f_sand)

      ! Volumetric soil water concentration at saturation point theta_s
      ! JULES_SOIL_PROPS units: m^3 water per m^3 soil
      sm_sat =  0.505 - 0.037 * f_clay - 0.142 * f_sand
        
      ! Volumetric soil moisture concentration at wilting point theta_w
      ! This is calculated assuming a matrix water potential (psi) of 1.5MPa
      ! 1.0 MPa = 1.0E6 kg m^-1 s^-2
      ! JULES_SOIL_PROPS units: m^3 water per m^3 soil
      sm_wilt = brooks_and_corey_equation(sm_sat, sathh, b,  1.5E6)
      
      ! Volumetric soil moisture concentration at critical point theta_c
      ! This is calculated assuming a matrix water potential (psi)
      ! of 0.033MPa.
      ! JULES_SOIL_PROPS units: m^3 water per m^3 soil
      sm_crit = brooks_and_corey_equation(sm_sat, sathh, b, 0.033E6)


      ! Calculate thermal conductivity of non-organic soil
      ! using [AncilDoc_SoilThermal] Method 1

      lambda_air  = 0.025  ! W m^-1 K^-1
      lambda_clay = 1.16   ! W m^-1 K^-1
      lambda_sand = 1.57   ! W m^-1 K^-1
      lambda_silt = 1.57   ! W m^-1 K^-1

      aa = (1.0 - sm_sat) * f_clay
      bb = (1.0 - sm_sat) * f_sand
      cc = (1.0 - sm_sat) * f_silt

      hcon = (lambda_air ** sm_sat * lambda_clay ** aa * lambda_sand ** bb * lambda_silt ** cc)
        
      !  Calculate heat capacity of non-organic soil
      !  using [AncilDoc_SoilThermal] Method 1 '''

      c_clay = 2.373E6  ! J m^-3 K^-1
      c_sand = 2.133E6  ! J m^-3 K^-1
      c_silt = 2.133E6  ! J m^-3 K^-1

      hcap = ((1.0 - sm_sat) * (f_clay * c_clay + f_sand * c_sand + f_silt * c_silt))
    else
      b        = -9999.0
      sathh    = -9999.0
      satcon   = -9999.0
      sm_sat   = -9999.0
      sm_wilt  = -9999.0
      sm_crit  = -9999.0
      hcap     = -9999.0
      hcon     = -9999.0
    endif
  end subroutine calc_soil_params_id
  
  function brooks_and_corey_equation(sm_sat, sathh, b,  psi) result(theta)
    implicit none
    real    :: sm_sat, sathh, b, psi 
    real    :: theta
    real    :: rho_w, g, h

    rho_w = 999.97  ! density of water kg m-3 (at 4 degC)
    g     = 9.81    ! acceleration due to gravity in m s-2

    ! h is the soil water pressure head
    h = psi / (rho_w * g)

    theta = sm_sat * ((sathh / h) ** (1.0 / b))
  end function brooks_and_corey_equation
    

end module JULES50_parmsMod

