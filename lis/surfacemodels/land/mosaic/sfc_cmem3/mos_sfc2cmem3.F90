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
!BOP
! !ROUTINE: mos_sfc2cmem3
! \label{mos_sfc2cmem3}
!
! !REVISION HISTORY:
!  37 Mar 2009: Sujay Kumar; Initial Code
!  2  May 2011: Yudong Tian; Modified for Mos 3.2 
!
! !INTERFACE:
subroutine mos_sfc2cmem3(n, sfcState)
! !USES:      
  use ESMF
  use LIS_coreMod,   only : LIS_rc, LIS_surface
  use LIS_logMod,    only : LIS_verify
  use LIS_soilsMod, only : LIS_soils
  use LIS_historyMod, only : LIS_patch2tile
  use LIS_constantsMod,  only : LIS_CONST_RHOFW

  use mos_lsmMod

  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: n
  type(ESMF_State)    :: sfcState
  
! FUNCTIONS

! 
! !DESCRIPTION: 
! This subroutine assigns the mos specific surface variables
! to CMEM3.
!
!EOP

  integer             :: t, status
  real                :: smc, wilt, por
  type(ESMF_Field)    :: windField, lcField, ltField, ltempField, stempField, &
       smcField, relsmcField, stcField, cmcField, vegFracField, &
       scField, snowhField, snodField, laiField
  real, pointer       :: wind_speed(:), land_coverage(:), land_temperature(:),  &
       soil_moisture_content(:), soil_temperature(:), canopy_water_content(:), &
       vegetation_fraction(:), snow_depth(:), snow_density(:), land_type(:), &
       snow_coverage(:), snow_temperature(:), leaf_area_index(:), relative_soil_moisture_content(:)
  real, allocatable    :: lsm_var(:)
  real, allocatable    :: snow_cov(:)
  real, allocatable    :: land_cov(:)

  allocate(lsm_var(LIS_rc%npatch(n,LIS_rc%lsm_index)))
  allocate(snow_cov(LIS_rc%npatch(n,LIS_rc%lsm_index)))
  allocate(land_cov(LIS_rc%npatch(n,LIS_rc%lsm_index)))

  call ESMF_StateGet(sfcState,"Wind Speed",windField,rc=status)
  call LIS_verify(status)
  call ESMF_FieldGet(windField,localDE=0,farrayPtr=wind_speed,rc=status)
  call LIS_verify(status)

  call ESMF_StateGet(sfcState,"Land Coverage",lcField,rc=status)
  call LIS_verify(status)
  call ESMF_FieldGet(lcField,localDE=0,farrayPtr=land_coverage,rc=status)
  call LIS_verify(status)

  call ESMF_StateGet(sfcState,"Land Type",ltField,rc=status)
  call LIS_verify(status)
  call ESMF_FieldGet(ltField,localDE=0,farrayPtr=land_type,rc=status)
  call LIS_verify(status)

  call ESMF_StateGet(sfcState,"Land Temperature",ltempField,rc=status)
  call LIS_verify(status)
  call ESMF_FieldGet(ltempField,localDE=0,farrayPtr=land_temperature,rc=status)
  call LIS_verify(status)

  call ESMF_StateGet(sfcState,"Snow Temperature",stempField,rc=status)
  call LIS_verify(status)
  call ESMF_FieldGet(stempField,localDE=0,farrayPtr=snow_temperature,rc=status)
  call LIS_verify(status)

  call ESMF_StateGet(sfcState,"Soil Moisture Content",smcField,rc=status)
  call LIS_verify(status)
  call ESMF_FieldGet(smcField,localDE=0,farrayPtr=soil_moisture_content,rc=status)
  call LIS_verify(status)

  call ESMF_StateGet(sfcState,"Soil Temperature",stcField,rc=status)
  call LIS_verify(status)
  call ESMF_FieldGet(stcField,localDE=0,farrayPtr=soil_temperature, rc=status)
  call LIS_verify(status)

  call ESMF_StateGet(sfcState,"Canopy Water Content",cmcField,rc=status)
  call LIS_verify(status)
  call ESMF_FieldGet(cmcField,localDE=0,farrayPtr=canopy_water_content,rc=status)
  call LIS_verify(status)

  call ESMF_StateGet(sfcState,"Vegetation Fraction",vegfracField,rc=status)
  call LIS_verify(status)
  call ESMF_FieldGet(vegfracField,localDE=0,farrayPtr=vegetation_fraction,rc=status)
  call LIS_verify(status)

  call ESMF_StateGet(sfcState,"Snow Coverage",scField,rc=status)
  call LIS_verify(status)
  call ESMF_FieldGet(scField,localDE=0,farrayPtr=snow_coverage,rc=status)
  call LIS_verify(status)

  call ESMF_StateGet(sfcState,"Snow Depth",snowhField,rc=status)
  call LIS_verify(status)
  call ESMF_FieldGet(snowhField,localDE=0,farrayPtr=snow_depth,rc=status)
  call LIS_verify(status)

  call ESMF_StateGet(sfcState,"Snow Density",snodField,rc=status)
  call LIS_verify(status)
  call ESMF_FieldGet(snodField,localDE=0,farrayPtr=snow_density,rc=status)
  call LIS_verify(status)

  call ESMF_StateGet(sfcState,"Leaf Area Index",laiField,rc=status)
  call LIS_verify(status)
  call ESMF_FieldGet(laiField,localDE=0,farrayPtr=leaf_area_index,rc=status)
  call LIS_verify(status)

  call ESMF_StateGet(sfcState,"Relative Soil Moisture Content",relsmcField,rc=status)
  call LIS_verify(status)
  call ESMF_FieldGet(relsmcField,localDE=0,farrayPtr=relative_soil_moisture_content,rc=status)
  call LIS_verify(status)

  !REWRITE THIS ROUTINE AS AN ESMF_STATE + FIELD STYLE. 
  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     ! obviously much work needs to be done here. 
     !gross type of surface determined by coverage
     ! assume only land points
     ! Surface type independent data
     lsm_var(t) = sqrt(mos_struc(n)%mos(t)%uwind*&
          mos_struc(n)%mos(t)%uwind+&
          mos_struc(n)%mos(t)%vwind*mos_struc(n)%mos(t)%vwind)
  enddo
  !put in wind direction using math on u and v; deg. E. of N.
  call LIS_patch2tile(n,LIS_rc%lsm_index,wind_speed,lsm_var)
     ! Snow surface type data
     ! Snow depth are m in CMEM3 and in mos
  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     lsm_var(t) = mos_struc(n)%mos(t)%snow*10.0/1000.0 
  enddo
  call LIS_patch2tile(n,LIS_rc%lsm_index,snow_depth,lsm_var)

  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     if(mos_struc(n)%mos(t)%snow.gt.0)  then 
        !Units are g/cm^3 in CMEM3.  

        lsm_var(t) =0.10

        ! intuitively: for 1kg water equivalent, if snowh high, density low

     else
        lsm_var(t) = 0.0
     endif
  enddo
  call LIS_patch2tile(n,LIS_rc%lsm_index,snow_density,lsm_var)
     ! Snow logic: 
     !   if snow exceeds some threshold, then assign 100% snow, else 100% land
     !   CHECK THRESHOLD VALIDITY.  SET EQUAL TO WEIZHONGS.  0.1 MM (0.0001 m).
  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     if (mos_struc(n)%mos(t)%snow  .gt. 0.0001) then  
        snow_cov(t)=1.0
        land_cov(t)=0.0
     else
        snow_cov(t)=0.0
        land_cov(t) = 1.0
     endif
  enddo
  call LIS_patch2tile(n,LIS_rc%lsm_index,snow_coverage,snow_cov)
  call LIS_patch2tile(n,LIS_rc%lsm_index,land_coverage,land_cov)

  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     lsm_var(t)     =  LIS_surface(n,LIS_rc%lsm_index)%tile(t)%vegt
  enddo
  call LIS_patch2tile(n,LIS_rc%lsm_index,land_type, lsm_var)
  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     lsm_var(t) = mos_struc(n)%mos(t)%ct
  enddo
  call LIS_patch2tile(n,LIS_rc%lsm_index,land_temperature, lsm_var)
  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     lsm_var(t) = mos_struc(n)%mos(t)%ct
  enddo
  call LIS_patch2tile(n,LIS_rc%lsm_index,snow_temperature, lsm_var)
  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     lsm_var(t) = mos_struc(n)%mos(t)%water1/&
          (mos_struc(n)%mos(t)%soilp(4)*LIS_CONST_RHOFW)
  enddo
  call LIS_patch2tile(n,LIS_rc%lsm_index,soil_moisture_content, lsm_var)
  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     lsm_var(t) = mos_struc(n)%mos(t)%soT
  enddo
  call LIS_patch2tile(n,LIS_rc%lsm_index,soil_temperature, lsm_var)
  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     lsm_var(t) = MOS_STRUC(N)%MOS(t)%VEGIP(1)
  enddo
  call LIS_patch2tile(n,LIS_rc%lsm_index,vegetation_fraction, lsm_var)
  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     ! --------------------------------------------------------
     ! Canopy water content (mm): only includes canopy interception (cmc). 
     ! VWC will be calculated in CMEM module with lai  
     ! --------------------------------------------------------
     lsm_var(t) = mos_struc(n)%mos(t)%ics*LIS_CONST_RHOFW 
  enddo
  call LIS_patch2tile(n,LIS_rc%lsm_index,canopy_water_content, lsm_var)
  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     lsm_var(t) = MOS_STRUC(N)%MOS(t)%VEGIP(2)
  enddo
  call LIS_patch2tile(n,LIS_rc%lsm_index,leaf_area_index, lsm_var)
  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     smc = mos_struc(n)%mos(t)%water1/&
          (mos_struc(n)%mos(t)%soilp(4)*LIS_CONST_RHOFW) 
     wilt = ((mos_struc(n)%mos(t)%vegp(11)/&
          mos_struc(n)%mos(t)%soilp(2)) ** (1.0 / &
          (-mos_struc(n)%mos(t)%soilp(1))))
     por  =  lis_soils(n)%porosity(&
          LIS_surface(n,LIS_rc%lsm_index)%tile(t)%col,&
          LIS_surface(n,LIS_rc%lsm_index)%tile(t)%row,1)
     lsm_var(t) =mos_struc(n)%mos(t)%sowet(1)
     if ( lsm_var(t) > 1.0 ) then
        lsm_var(t)  = 1.0
     endif
     if ( lsm_var(t)  < 0.01 ) then
        lsm_var(t)  = 0.01
     endif
  enddo
  call LIS_patch2tile(n,LIS_rc%lsm_index,relative_soil_moisture_content, lsm_var)

  deallocate(lsm_var)
  deallocate(snow_cov)
  deallocate(land_cov)

end subroutine mos_sfc2cmem3


