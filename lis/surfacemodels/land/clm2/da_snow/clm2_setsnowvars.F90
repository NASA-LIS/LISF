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
! !ROUTINE: clm2_setsnowvars
!  \label{clm2_setsnowvars}
!
! !REVISION HISTORY:
! 27Feb2005: Sujay Kumar; Initial Specification
! 25Jun2006: Sujay Kumar: Updated for the ESMF design
! 06Oct2008: Gabrielle De Lannoy: update prognostic 
!            and directly derived diagnostic vars for sca/f EnKFiltering
! 04Nov2008: Gabrielle De Lannoy: only for total snow water
!
! !INTERFACE:
subroutine clm2_setsnowvars(n, LSM_State)
! !USES:
  use ESMF
  use LIS_coreMod, only : LIS_rc
  use LIS_logMod,  only  : LIS_verify
  use clm2_lsmMod
  use clm2_varpar,    only : nlevsno
  use clm2_varcon,    only : denh2o, tfrz, bdsno

  implicit none
! !ARGUMENTS: 
  integer, intent(in)    :: n
  type(ESMF_State)       :: LSM_State
!
! !DESCRIPTION:
!  
!  This routine assigns the soil moisture prognostic variables to clm's
!  model space. 
! 
!EOP

  type(ESMF_Field)       :: h2osnoField
  type(ESMF_Field)       :: snowdpField

  real, pointer          :: h2osno(:)  ! clm2.0 total snow water (kg m^-2)
  real, pointer          :: snowdp(:)  ! clm2.0 total snow depth

  integer                :: t,j
  integer                :: status
  real                   :: incr_h2osno, incr_snowdp
  logical                :: snow_check

! -----------------------------------------------------------------------------

  call ESMF_StateGet(LSM_State,"Total Snow Water",h2osnoField,rc=status)
  call LIS_verify(status)
  call ESMF_StateGet(LSM_State,"Total Snow Depth",snowdpField,rc=status)
  call LIS_verify(status)

  call ESMF_FieldGet(h2osnoField,localDE=0,farrayPtr=h2osno,rc=status)
  call LIS_verify(status)
  call ESMF_FieldGet(snowdpField,localDE=0,farrayPtr=snowdp,rc=status)
  call LIS_verify(status)

  do t=1,LIS_rc%ntiles(n)
     incr_h2osno = clm2_struc(n)%clm(t)%h2osno
     incr_snowdp = clm2_struc(n)%clm(t)%snowdp

     if(clm2_struc(n)%clm(t)%do_capsnow) then 
        clm2_struc(n)%clm(t)%h2osno = min(h2osno(t),clm2_struc(n)%clm(t)%h2osno)
        clm2_struc(n)%clm(t)%snowdp = min(snowdp(t),clm2_struc(n)%clm(t)%snowdp)
     else
        clm2_struc(n)%clm(t)%h2osno = h2osno(t)
        clm2_struc(n)%clm(t)%snowdp = snowdp(t)
     endif

   ! Re-initialize depths and water content of indiv layer, only ice (as in Hydrology1)
   !  simply use a fraction to update all sub-state-variables:

     if(clm2_struc(n)%clm(t)%snl < 0 ) then 
        
        incr_h2osno = clm2_struc(n)%clm(t)%h2osno / incr_h2osno
        incr_snowdp = clm2_struc(n)%clm(t)%snowdp / incr_snowdp

        do j = clm2_struc(n)%clm(t)%snl+1, 0
           clm2_struc(n)%clm(t)%h2osoi_ice(j) = clm2_struc(n)%clm(t)%h2osoi_ice(j)*incr_h2osno
           clm2_struc(n)%clm(t)%dz(j) = clm2_struc(n)%clm(t)%dz(j)*incr_snowdp
        enddo
        snow_check = .false. 
        do j=1,clm2_struc(n)%clm(t)%snl+1
           snow_check = snow_check.or.(clm2_struc(n)%clm(t)%dz(j).ne.0)
        enddo

        if(snow_check) then 
           call SnowCompaction (clm2_struc(n)%clm(t))
           call CombineSnowLayers (clm2_struc(n)%clm(t))
           call DivideSnowLayers (clm2_struc(n)%clm(t))
           
            if (clm2_struc(n)%clm(t)%snl > -nlevsno) then
              clm2_struc(n)%clm(t)%snowage = 0.
              do j = -nlevsno+1, clm2_struc(n)%clm(t)%snl
                 clm2_struc(n)%clm(t)%h2osoi_ice(j) = 0.
                 clm2_struc(n)%clm(t)%h2osoi_liq(j) = 0.
                 clm2_struc(n)%clm(t)%t_soisno(j)   = 0.
                 clm2_struc(n)%clm(t)%dz(j)         = 0.
                 clm2_struc(n)%clm(t)%z(j)          = 0.
                 clm2_struc(n)%clm(t)%zi(j-1)       = 0.
              enddo
           endif
        
        elseif (clm2_struc(n)%clm(t)%snowdp > 0.01 .and. clm2_struc(n)%clm(t)%h2osno > 0.01*bdsno) then
           clm2_struc(n)%clm(t)%snl            = -1
           clm2_struc(n)%clm(t)%dz(0)          = clm2_struc(n)%clm(t)%snowdp              ! meter
           clm2_struc(n)%clm(t)%z(0)           = -0.5*clm2_struc(n)%clm(t)%dz(0)
           clm2_struc(n)%clm(t)%zi(-1)         = -clm2_struc(n)%clm(t)%dz(0)
           clm2_struc(n)%clm(t)%snowage        = 0.                                      ! snow age
           clm2_struc(n)%clm(t)%t_soisno (0)   = min(tfrz, clm2_struc(n)%clm(t)%forc_t)   ! K
           clm2_struc(n)%clm(t)%h2osoi_ice(0)  = clm2_struc(n)%clm(t)%h2osno              ! kg m-2
           clm2_struc(n)%clm(t)%h2osoi_liq(0)  = 0.                                      ! kg m-2
           clm2_struc(n)%clm(t)%frac_iceold(0) = 1.        
        else
           clm2_struc(n)%clm(t)%snowage = 0.
           clm2_struc(n)%clm(t)%snl = 0
        endif
     endif
  enddo
end subroutine clm2_setsnowvars

