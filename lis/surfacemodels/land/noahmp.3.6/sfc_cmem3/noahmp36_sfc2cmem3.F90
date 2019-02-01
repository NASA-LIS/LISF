!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LIS_misc.h"
!BOP
! !ROUTINE: noahmp36_sfc2cmem3
! \label{noahmp36_sfc2cmem3}
!
! !REVISION HISTORY:
!  37 Mar 2009: Sujay Kumar; Initial Code
!  2  May 2011: Yudong Tian; Modified for Noah 3.2 
!  22 Jun 2018: Peter Shellito; Modified for NoahMP 3.6. 
!
! !INTERFACE:
subroutine noahmp36_sfc2cmem3(n, sfcState)
! !USES:      
  use ESMF
  use LIS_coreMod,   only : LIS_rc
  use LIS_logMod,    only : LIS_verify
  use LIS_historyMod, only : LIS_patch2tile
  use LIS_constantsMod,  only : LIS_CONST_RHOFW

  use NoahMP36_lsmMod

  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: n
  type(ESMF_State)    :: sfcState
  
! FUNCTIONS

! 
! !DESCRIPTION: 
! This subroutine assigns the noahMP36 specific surface variables
! to CMEM3.
!
!EOP

! Initialize variables
  integer             :: t, status
! These vars are of type ESMF_Field
! The next 3 lines contain the fields but without relsmcField
!  type(ESMF_Field)    :: windField, lcField, ltField, ltempField, stempField, &
!       smcField, stcField, cmcField, vegFracField, &
!       scField, snowhField, snodField, laiField
  type(ESMF_Field)    :: windField, lcField, ltField, ltempField, stempField, &
       smcField, relsmcField, stcField, cmcField, vegFracField, &
       scField, snowhField, snodField, laiField
! These vars are pointers
  real, pointer       :: wind_speed(:), land_coverage(:), land_temperature(:), &
       soil_moisture_content(:), soil_temperature(:), canopy_water_content(:), &
       vegetation_fraction(:), snow_depth(:), snow_density(:), land_type(:), &
       snow_coverage(:), snow_temperature(:), leaf_area_index(:), &
       relative_soil_moisture_content(:)
!       snow_coverage(:), snow_temperature(:), leaf_area_index(:)
! the line above doesn't have relative_soil_moisture_content
  real, allocatable       :: lsm_var(:)
  real, allocatable       :: snow_cov(:)
  real, allocatable       :: land_cov(:)

! And here we assign/determine the actual bounds for the array
  allocate(lsm_var(LIS_rc%npatch(n,LIS_rc%lsm_index)))
  allocate(snow_cov(LIS_rc%npatch(n,LIS_rc%lsm_index)))
  allocate(land_cov(LIS_rc%npatch(n,LIS_rc%lsm_index)))

! I don't have to do anything with these lines. 
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
  call ESMF_FieldGet(smcField,localDE=0,farrayPtr=soil_moisture_content,&
       rc=status)
  call LIS_verify(status)

  call ESMF_StateGet(sfcState,"Soil Temperature",stcField,rc=status)
  call LIS_verify(status)
  call ESMF_FieldGet(stcField,localDE=0,farrayPtr=soil_temperature, rc=status)
  call LIS_verify(status)

  call ESMF_StateGet(sfcState,"Canopy Water Content",cmcField,rc=status)
  call LIS_verify(status)
  call ESMF_FieldGet(cmcField,localDE=0,farrayPtr=canopy_water_content,&
       rc=status)
  call LIS_verify(status)

  call ESMF_StateGet(sfcState,"Vegetation Fraction",vegfracField,rc=status)
  call LIS_verify(status)
  call ESMF_FieldGet(vegfracField,localDE=0,farrayPtr=vegetation_fraction,&
       rc=status)
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

  call ESMF_StateGet(sfcState,"Relative Soil Moisture Content",relsmcField,&
       rc=status)
  call LIS_verify(status)
  call ESMF_FieldGet(relsmcField,localDE=0,&
       farrayPtr=relative_soil_moisture_content,rc=status)
  call LIS_verify(status)

! In this section, we loop through the "patches" to accumulate the variable for each patch, then the routine LIS_patch2tile accumulates that patch over each tile and saves it as the variable named in the second-to-last parameter there.
! The LIS_patch2tile is found in LIS_historyMod.F90 in /core
! public :: LIS_patch2tile ! convert a patchspace variable to tilespace
! subroutine LIS_patch2tile(n,m,tvar,pvar):
! n is the index of the domain or nest
! m is unknown, probably another index
! tvar is the variable dimensioned in the tile space
! pvar is the variable dimensioned in the patch space

  !REWRITE THIS ROUTINE AS AN ESMF_STATE + FIELD STYLE. 
  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     ! obviously much work needs to be done here. 
     !gross type of surface determined by coverage
     ! assume only land points
     ! Surface type independent data
     lsm_var(t) = sqrt(NOAHMP36_struc(n)%noahmp36(t)%wind_e*& 
          NOAHMP36_struc(n)%noahmp36(t)%wind_e+& 
          NOAHMP36_struc(n)%noahmp36(t)%wind_n*&
          NOAHMP36_struc(n)%noahmp36(t)%wind_n)
  enddo
  !put in wind direction using math on u and v; deg. E. of N.
  call LIS_patch2tile(n,LIS_rc%lsm_index,wind_speed,lsm_var)
     ! Snow surface type data
     ! Snow depth are m in CMEM3 and in noah
  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     lsm_var(t) = NOAHMP36_struc(n)%noahmp36(t)%snowh 
  enddo
  call LIS_patch2tile(n,LIS_rc%lsm_index,snow_depth,lsm_var)
  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     if(NOAHMP36_struc(n)%noahmp36(t)%sneqv.gt.0)  then 
        ! CMEM3 requires snow density in units of g/cm^3 in CMEM3.
        ! snowh is in 'm snow', sneqv is in 'mm water'
        ! divide sneqv by 1000 to convert it to m water.
        !Density of water: 1 g/cm^3
        !By definition of sneqv:
        ! mass of snow per unit area     = mass of snow water equivalent per
        !                                     unit area
        ! snowh*snowdensity = sneqv*water density
        ! snowdensity = (sneqv/snowh)*water density

        lsm_var(t) = NOAHMP36_struc(n)%noahmp36(t)%sneqv &
             / NOAHMP36_struc(n)%noahmp36(t)%snowh / 1000.0

     else
        lsm_var(t) = 0.0
     endif
  enddo
  call LIS_patch2tile(n,LIS_rc%lsm_index,snow_density,lsm_var)
     ! Snow logic: 
     !   if snow exceeds some threshold, then assign 100% snow, else 100% land
     !   CHECK THRESHOLD VALIDITY.  SET EQUAL TO WEIZHONGS.  0.1 MM (0.0001 m).
  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     if (NOAHMP36_struc(n)%noahmp36(t)%snowh  .gt. 0.0001) then  
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
     lsm_var(t)     = NOAHMP36_struc(n)%noahmp36(t)%vegetype 
  enddo
  call LIS_patch2tile(n,LIS_rc%lsm_index,land_type, lsm_var)
  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     lsm_var(t) = NOAHMP36_struc(n)%noahmp36(t)%trad
  enddo
  call LIS_patch2tile(n,LIS_rc%lsm_index,land_temperature, lsm_var)
  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     lsm_var(t) = NOAHMP36_struc(n)%noahmp36(t)%trad
  enddo
  call LIS_patch2tile(n,LIS_rc%lsm_index,snow_temperature, lsm_var)
  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     lsm_var(t) = NOAHMP36_struc(n)%noahmp36(t)%smc(1)
  enddo
  call LIS_patch2tile(n,LIS_rc%lsm_index,soil_moisture_content, lsm_var)
  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     lsm_var(t) = NOAHMP36_struc(n)%noahmp36(t)%sstc(NOAHMP36_struc(n)%nsnow+1)
  enddo
  call LIS_patch2tile(n,LIS_rc%lsm_index,soil_temperature, lsm_var)
  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     lsm_var(t) = NOAHMP36_struc(n)%noahmp36(t)%fveg
  enddo
  call LIS_patch2tile(n,LIS_rc%lsm_index,vegetation_fraction, lsm_var)
  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     ! --------------------------------------------------------
     ! Canopy water content (mm): only includes canopy interception (cmc). 
     ! VWC will be calculated in CMEM module with lai  
     ! --------------------------------------------------------
     lsm_var(t) = NOAHMP36_struc(n)%noahmp36(t)%canliq
  enddo
  call LIS_patch2tile(n,LIS_rc%lsm_index,canopy_water_content, lsm_var)
  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     lsm_var(t) = NOAHMP36_struc(n)%noahmp36(t)%lai
  enddo
  call LIS_patch2tile(n,LIS_rc%lsm_index,leaf_area_index, lsm_var)
  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
           lsm_var(t) = &  
                (NOAHMP36_struc(n)%noahmp36(t)%smc(1) - &
                NOAHMP36_struc(n)%noahmp36(t)%smcwlt) / &
                (NOAHMP36_struc(n)%noahmp36(t)%smcmax - &
                NOAHMP36_struc(n)%noahmp36(t)%smcwlt)
           if ( lsm_var(t) > 1.0 ) then
              lsm_var(t)  = 1.0
           endif
           if ( lsm_var(t)  < 0.01 ) then
              lsm_var(t)  = 0.01
           endif
  enddo
  call LIS_patch2tile(n,LIS_rc%lsm_index,relative_soil_moisture_content, &
       lsm_var)

  deallocate(lsm_var)
  deallocate(snow_cov)
  deallocate(land_cov)

end subroutine noahmp36_sfc2cmem3
