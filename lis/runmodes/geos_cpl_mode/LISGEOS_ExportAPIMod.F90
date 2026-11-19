!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.8
!
! Copyright (c) 2026 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------

module LISGEOS_ExportAPIMod
!
!BOP
! !MODULE: LISGEOS_ExportAPIMod
!
! !DESCRIPTION:
!  This module provides a public API for the GEOS MAPL cap to extract
!  explicitly defined LIS export variables for coupling. It grabs
!  data directly from the pre-populated LISGEOS_export structure.
!
! !USES:
  use LISGEOSGridCompMod, only : LISGEOS_export
!EOP

  implicit none
  private

  public :: LISGEOS_GetExportVar

contains

!BOP
! !ROUTINE: LISGEOS_GetExportVar
!
! !DESCRIPTION:
!  Retrieves a pre-calculated 2D LIS export grid by variable name.
!  Returns status = 0 on success, or 1 if the variable is not found.
!
! !INTERFACE:
  subroutine LISGEOS_GetExportVar(n, var_name, gridvar, status)
! !ARGUMENTS:
    implicit none
    integer,          intent(in)    :: n           ! LIS Nest index
    character(len=*), intent(in)    :: var_name    ! Name of variable requested by GEOS
    real,             intent(inout) :: gridvar(:,:)! 2D Grid array pointer provided by GEOS
    integer,          intent(out)   :: status      ! Return status (0=Success, 1=Error)
!EOP

    ! Initialize success
    status = 0

    ! Map the requested variable name to the explicit LISGEOS_export 2D array
    select case (trim(var_name))
       
       ! ==========================================
       ! 1. FLUXES, ENERGY, & BASIC SURFACE STATES
       ! ==========================================
       case ('Qh', 'QH', 'qh')
          gridvar(:,:) = LISGEOS_export(n)%qh(:,:)
       case ('Qle', 'QLE', 'qle')
          gridvar(:,:) = LISGEOS_export(n)%qle(:,:)
       case ('AvgSurfT', 'AVGSURFT', 'avgsurft')
          gridvar(:,:) = LISGEOS_export(n)%avgsurft(:,:)
       case ('Emiss', 'EMISS', 'emiss')
          gridvar(:,:) = LISGEOS_export(n)%emiss(:,:)
       case ('Albedo', 'ALBEDO', 'albedo')
          gridvar(:,:) = LISGEOS_export(n)%albedo(:,:)
       case ('Qs', 'QS', 'qs')
          gridvar(:,:) = LISGEOS_export(n)%qs(:,:)
       case ('Qsb', 'QSB', 'qsb')
          gridvar(:,:) = LISGEOS_export(n)%qsb(:,:)
       case ('RootMoist', 'ROOTMOIST', 'rootmoist')
          gridvar(:,:) = LISGEOS_export(n)%rootmoist(:,:)
       case ('Soilm', 'SOILM', 'soilm')
          gridvar(:,:) = LISGEOS_export(n)%soilm(:,:)
       case ('TotalRunoff', 'TOTALRUNOFF', 'totalrunoff')
          gridvar(:,:) = LISGEOS_export(n)%totalrunoff(:,:)

#ifdef WRF_HYDRO
       case ('Infxsrt', 'INFXSRT', 'infxsrt')
          gridvar(:,:) = LISGEOS_export(n)%infxsrt(:,:)
#endif
          
       ! ==========================================
       ! 2. SOIL LAYER FIELDS
       ! ==========================================
       ! Soil Moisture Content
       case ('Smc1', 'SMC1', 'smc1')
          gridvar(:,:) = LISGEOS_export(n)%smc1(:,:)
       case ('Smc2', 'SMC2', 'smc2')
          gridvar(:,:) = LISGEOS_export(n)%smc2(:,:)
       case ('Smc3', 'SMC3', 'smc3')
          gridvar(:,:) = LISGEOS_export(n)%smc3(:,:)
       case ('Smc4', 'SMC4', 'smc4')
          gridvar(:,:) = LISGEOS_export(n)%smc4(:,:)
       ! Volumetric Liquid Soil Moisture
       case ('Sh2o1', 'SH2O1', 'sh2o1')
          gridvar(:,:) = LISGEOS_export(n)%sh2o1(:,:)
       case ('Sh2o2', 'SH2O2', 'sh2o2')
          gridvar(:,:) = LISGEOS_export(n)%sh2o2(:,:)
       case ('Sh2o3', 'SH2O3', 'sh2o3')
          gridvar(:,:) = LISGEOS_export(n)%sh2o3(:,:)
       case ('Sh2o4', 'SH2O4', 'sh2o4')
          gridvar(:,:) = LISGEOS_export(n)%sh2o4(:,:)
       ! Soil Temperature
       case ('Stc1', 'STC1', 'stc1')
          gridvar(:,:) = LISGEOS_export(n)%stc1(:,:)
       case ('Stc2', 'STC2', 'stc2')
          gridvar(:,:) = LISGEOS_export(n)%stc2(:,:)
       case ('Stc3', 'STC3', 'stc3')
          gridvar(:,:) = LISGEOS_export(n)%stc3(:,:)
       case ('Stc4', 'STC4', 'stc4')
          gridvar(:,:) = LISGEOS_export(n)%stc4(:,:)

       ! ==========================================
       ! 3. SNOW FIELDS
       ! ==========================================
       case ('Snow', 'SNOW', 'snow')
          gridvar(:,:) = LISGEOS_export(n)%snow(:,:)
       case ('Snowh', 'SNOWH', 'snowh')
          gridvar(:,:) = LISGEOS_export(n)%snowh(:,:)
       case ('Snocvr', 'SNOCVR', 'snocvr')
          gridvar(:,:) = LISGEOS_export(n)%snocvr(:,:)
       case ('SnowDepth_SCA', 'SNOWDEPTH_SCA', 'snowdepth_sca')
          gridvar(:,:) = LISGEOS_export(n)%snowdepth_sca(:,:)

       ! ==========================================
       ! 4. AERODYNAMIC & SURFACE EXCHANGE
       ! ==========================================
       case ('Rhos', 'RHOS', 'rhos')
          gridvar(:,:) = LISGEOS_export(n)%rhos(:,:)
       case ('Ch_f', 'CH_F', 'ch_f')
          gridvar(:,:) = LISGEOS_export(n)%ch_f(:,:)
       case ('Cm_f', 'CM_F', 'cm_f')
          gridvar(:,:) = LISGEOS_export(n)%cm_f(:,:)
       case ('AvgSurfQ', 'AVGSURFQ', 'avgsurfq')
          gridvar(:,:) = LISGEOS_export(n)%avgsurfq(:,:)
       case ('Z0', 'z0')
          gridvar(:,:) = LISGEOS_export(n)%z0(:,:)
       case ('Zpd', 'ZPD', 'zpd')
          gridvar(:,:) = LISGEOS_export(n)%zpd(:,:)
       case ('UstarVeg', 'USTARVEG', 'ustarveg')
          gridvar(:,:) = LISGEOS_export(n)%ustarveg(:,:)
       case ('UstarBare', 'USTARBARE', 'ustarbare')
          gridvar(:,:) = LISGEOS_export(n)%ustarbare(:,:)
       case ('Taux', 'TAUX', 'taux')
          gridvar(:,:) = LISGEOS_export(n)%taux(:,:)
       case ('Tauy', 'TAUY', 'tauy')
          gridvar(:,:) = LISGEOS_export(n)%tauy(:,:)

       ! ==========================================
       ! 5. ADVANCED ALBEDO FIELDS
       ! ==========================================
       case ('AlbedoVisDir', 'ALBEDOVISDIR', 'albedovisdir')
          gridvar(:,:) = LISGEOS_export(n)%albedovisdir(:,:)
       case ('AlbedoNirDir', 'ALBEDONIRDIR', 'albedonirdir')
          gridvar(:,:) = LISGEOS_export(n)%albedonirdir(:,:)
       case ('AlbedoVisDif', 'ALBEDOVISDIF', 'albedovisdif')
          gridvar(:,:) = LISGEOS_export(n)%albedovisdif(:,:)
       case ('AlbedoNirDif', 'ALBEDONIRDIF', 'albedonirdif')
          gridvar(:,:) = LISGEOS_export(n)%albedonirdif(:,:)

       ! ==========================================
       ! 6. VEGETATION, CARBON, & MASS POOLS
       ! ==========================================
       case ('Rtmass', 'RTMASS', 'rtmass')
          gridvar(:,:) = LISGEOS_export(n)%rtmass(:,:)
       case ('Stmass', 'STMASS', 'stmass')
          gridvar(:,:) = LISGEOS_export(n)%stmass(:,:)
       case ('Wood', 'WOOD', 'wood')
          gridvar(:,:) = LISGEOS_export(n)%wood(:,:)
       case ('Stblcp', 'STBLCP', 'stblcp')
          gridvar(:,:) = LISGEOS_export(n)%stblcp(:,:)
       case ('Fastcp', 'FASTCP', 'fastcp')
          gridvar(:,:) = LISGEOS_export(n)%fastcp(:,:)
       case ('Lai', 'LAI', 'lai')
          gridvar(:,:) = LISGEOS_export(n)%lai(:,:)
       case ('Sai', 'SAI', 'sai')
          gridvar(:,:) = LISGEOS_export(n)%sai(:,:)
       case ('Vai', 'VAI', 'vai')
          gridvar(:,:) = LISGEOS_export(n)%vai(:,:)
       case ('Evap', 'EVAP', 'evap')
          gridvar(:,:) = LISGEOS_export(n)%evap(:,:)
       case ('Greenness', 'GREENNESS', 'greenness')
          gridvar(:,:) = LISGEOS_export(n)%greenness(:,:)
       case ('Landcover', 'LANDCOVER', 'landcover')
          gridvar(:,:) = LISGEOS_export(n)%landcover(:,:)
       case ('Nee', 'NEE', 'nee')
          gridvar(:,:) = LISGEOS_export(n)%nee(:,:)
       case ('CanopInt', 'CANOPINT', 'canopint')
          gridvar(:,:) = LISGEOS_export(n)%canopint(:,:)
       case ('Tair_t2m', 'TAIR_T2M', 'tair_t2m')
          gridvar(:,:) = LISGEOS_export(n)%tair_t2m(:,:)

       ! ==========================================
       ! ERROR HANDLING
       ! ==========================================
       case default
          ! The variable name requested by GEOS is not plumbed in this dictionary
          status = 1
    end select

  end subroutine LISGEOS_GetExportVar

end module LISGEOS_ExportAPIMod