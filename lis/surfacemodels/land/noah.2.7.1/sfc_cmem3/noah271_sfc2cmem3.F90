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
! !ROUTINE: noah271_sfc2cmem3
! \label{noah271_sfc2cmem3}
!
! !REVISION HISTORY:
!  37 Mar 2009: Sujay Kumar; Initial Code
!
! !INTERFACE:
subroutine noah271_sfc2cmem3(n, sfcState)
! !USES:      
  use ESMF
  use LIS_coreMod,   only : LIS_rc
  use LIS_logMod,    only : LIS_verify
  use LIS_constantsMod,  only : LIS_CONST_RHOFW

  use noah271_lsmMod

  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: n
  type(ESMF_State)    :: sfcState
  
! FUNCTIONS

! 
! !DESCRIPTION: 
! This subroutine assigns the noah271 specific surface variables
! to CMEM3. 
!
!EOP

  integer             :: t, status
  type(ESMF_Field)    :: windField, lcField, ltField, ltempField, stempField, &
       smcField, stcField, cmcField, vegFracField, &
       scField, snowhField, snodField
  real, allocatable       :: wind_speed(:), land_coverage(:), land_temperature(:),  &
       soil_moisture_content(:), soil_temperature(:), canopy_water_content(:), &
       vegetation_fraction(:), snow_depth(:), snow_density(:), land_type(:), &
       snow_coverage(:), snow_temperature(:)

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


  !REWRITE THIS ROUTINE AS AN ESMF_STATE + FIELD STYLE. 
  do t=1,LIS_rc%ntiles(n)
     ! obviously much work needs to be done here. 
     !gross type of surface determined by coverage
     ! assume only land points


     ! Surface type independent data
     wind_speed(t) = sqrt(noah271_struc(n)%noah(t)%uwind*noah271_struc(n)%noah(t)%uwind+&
          noah271_struc(n)%noah(t)%vwind*noah271_struc(n)%noah(t)%vwind)
     !put in wind direction using math on u and v; deg. E. of N.

     ! Snow surface type data

     !Units are m in CMEM3 and in NOAH: 
     Snow_Depth(t) = noah271_struc(n)%noah(t)%snowh

        ! Snow Density: 
        !Units are g/cm^3 in CMEM3.  snowh is in 'm snow', sneqv in 'm water'
        !Density H20 water: 1000 kg/m^3, or 1000 kg/m^3*1000g/kg = 10^6 g/m^3
        !By definition of sneqv:
        ! mass of snow per unit area     = mass of snow water equivalent per
        !                                     unit area
        ! 1*snowh*snowdensity = 1*sneqv*water density
        ! snowdensity = (sneqv/snowh)*water density

     if(noah271_struc(n)%noah(t)%sneqv.gt.0)  then 
        Snow_Density(t) =noah271_struc(n)%noah(t)%sneqv &
             / noah271_struc(n)%noah(t)%snowh

     else
        snow_density(t) = 0.0
     endif

     ! Snow logic: 
     !   if snow exceeds some threshold, then assign 100% snow, else 100% land
     !   CHECK THRESHOLD VALIDITY.  SET EQUAL TO WEIZHONGS.  0.1 MM (0.0001 m).
     if (Snow_Depth(t) .gt. 0.0001) then  
        Snow_Coverage(t)=1.0
        Land_Coverage(t)=0.0
     else
        Snow_Coverage(t)=0.0
        Land_Coverage(t) = 1.0
     endif
     Land_Type(t)     = noah271_struc(n)%noah(t)%vegt 
     Land_Temperature(t) = noah271_struc(n)%noah(t)%t1
     Snow_Temperature(t) = noah271_struc(n)%noah(t)%t1
     soil_moisture_content(t) = noah271_struc(n)%noah(t)%smc(1)
     soil_temperature(t) = noah271_struc(n)%noah(t)%stc(1)
     Canopy_Water_Content(t) = noah271_struc(n)%noah(t)%cmc*LIS_CONST_RHOFW
     Vegetation_Fraction(t) = noah271_struc(n)%noah(t)%vegip
  enddo
end subroutine noah271_sfc2cmem3


