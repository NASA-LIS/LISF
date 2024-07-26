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
! !ROUTINE: iniTimeVar
! \label{iniTimeVar}
! 
! !INTERFACE:
subroutine iniTimeVar (n, eccen, obliqr, lambm0 , mvelpp)
! !USES:
  use LIS_precisionMod
  use LIS_coreMod, only : LIS_rc
  use clm2_lsmMod
  use LIS_timeMgrMod, only : LIS_get_curr_calday
  use clm2_varpar,   only : nlevsno, nlevsoi
  use clm2_varcon  , only : bdsno, istice, istwet, &
       istsoil, denice, denh2o, tfrz, spval, doalb
  use clm2_shr_sys_mod , only : clm2_shr_sys_abort
  use LIS_logMod, only : LIS_logunit

  implicit none
! !ARGUMENTS: 
  integer , intent(in) :: n
  real(r8), intent(in) :: eccen    
  real(r8), intent(in) :: obliqr   
  real(r8), intent(in) :: lambm0   
  real(r8), intent(in) :: mvelpp   

! 
! !DESCRIPTION: 
! 
! Initialize the following time varying variables:
!
!   water      : h2osno, h2ocan, h2osoi$_-$liq, h2osoi$_-$ice, h2osoi$_-$vol \newline
!   snow       : snowdp, snowage, snl, dz, z, zi \newline
!   temperature: t$_-$soisno, t$_-$veg, t$_-$grnd \newline
!
! Note - h2osoi$_-$vol is needed by clm$_-$soilalb -this is not needed on 
! restart since it is computed before the soil albedo computation is 
! called
! 
! Note -  remaining variables are initialized by calls to ecosystem 
! dynamics and albedo subroutines. 
!
! 
! The arguments are: 
! \begin{description}
! \item[n]
!  index of the nest
! \item[eccen]
!  Earth's orbital eccentricity
! \item[obliqr]
!  Earth's obliquity in radians
! \item[lambm0]
!  Mean longitude of perihelion at the vernal equinox (radians)
! \item[mvelpp]
!  Earth's moving vernal equinox long. of perihelion + pi (radians)
! \end{description} 
!EOP

! --------------------------------------------------------------------

! ------------------------ local variables ---------------------------
  integer  :: i,k              !loop indices
  real(r8) :: calday                    !calendar day
! --------------------------------------------------------------------

! ----------------------------------------------------------------------
! Initialize water and temperature based on:
! ----------------------------------------------------------------------

  write (LIS_logunit,*) 'Setting initial data to non-spun up values'

! ========================================================================
! Set snow water 
! ========================================================================

! NOTE: h2ocan, h2osno, snowdp and snowage has valid values everywhere

  do k = 1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     clm2_struc(n)%clm(k)%h2ocan = 0.
     if (clm2_struc(n)%clm(k)%itypwat == istice) then
        clm2_struc(n)%clm(k)%h2osno = 1000.
     else
        clm2_struc(n)%clm(k)%h2osno = clm2_struc(n)%clm_iscv
     endif
     clm2_struc(n)%clm(k)%snowdp  = clm2_struc(n)%clm(k)%h2osno/bdsno
     clm2_struc(n)%clm(k)%snowage = 0.
  end do
     
! ========================================================================
! Set snow layer number, depth and thickiness 
! ========================================================================
  call snowdp2lev (n)
     
! ========================================================================
! Set snow/soil temperature
! ========================================================================
        
! NOTE: 
! t$_-$soisno only has valid values over non-lake
! t$_-$lake   only has valid values over lake
! t$_-$grnd has valid values over all land
! t$_-$veg  has valid values over all land

  do k =1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     clm2_struc(n)%clm(k)%t_soisno(-nlevsno+1:nlevsoi) = 0
     clm2_struc(n)%clm(k)%t_veg = clm2_struc(n)%clm_it
     if (.not. clm2_struc(n)%clm(k)%lakpoi) then  !not lake
        clm2_struc(n)%clm(k)%t_soisno(-nlevsno+1:0) = spval
        if (clm2_struc(n)%clm(k)%snl < 0) then    !snow layer temperatures
           do i = clm2_struc(n)%clm(k)%snl+1, 0     
              if (clm2_struc(n)%clm_it  < 273.16) then
                 clm2_struc(n)%clm(k)%t_soisno(i) = clm2_struc(n)%clm_it
              else
                 clm2_struc(n)%clm(k)%t_soisno(i) = 273.16 - 1.
              endif
           enddo
        endif
        do i = 1, nlevsoi
           if (clm2_struc(n)%clm(k)%itypwat == istice) then
              clm2_struc(n)%clm(k)%t_soisno(i) = clm2_struc(n)%clm_it
           else if (clm2_struc(n)%clm(k)%itypwat == istwet) then
              clm2_struc(n)%clm(k)%t_soisno(i) = clm2_struc(n)%clm_it
           else
              clm2_struc(n)%clm(k)%t_soisno(i) = clm2_struc(n)%clm_it
           endif
        enddo
        clm2_struc(n)%clm(k)%t_grnd = clm2_struc(n)%clm(k)%t_soisno(clm2_struc(n)%clm(k)%snl+1)
     else                           !lake
        clm2_struc(n)%clm(k)%t_grnd = clm2_struc(n)%clm_it
        clm2_struc(n)%clm(k)%t_lake = clm2_struc(n)%clm_it
     endif
  end do
  
! ========================================================================
! Set snow/soil ice and liquid mass
! ========================================================================
        
! volumetric water is set first and liquid content and ice lens are
! then obtained
! NOTE: h2osoi$_-$vol, h2osoi$_-$liq and h2osoi$_-$ice only have valid values 
! over soil
  do k = 1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     clm2_struc(n)%clm(k)%h2osoi_vol(         1:nlevsoi) = spval
     clm2_struc(n)%clm(k)%h2osoi_liq(-nlevsno+1:nlevsoi) = spval
     clm2_struc(n)%clm(k)%h2osoi_ice(-nlevsno+1:nlevsoi) = spval
     
     if (.not. clm2_struc(n)%clm(k)%lakpoi) then  !not lake
        if (clm2_struc(n)%clm(k)%snl < 0) then    !snow 
           do i = clm2_struc(n)%clm(k)%snl+1, 0
              clm2_struc(n)%clm(k)%h2osoi_ice(i) = clm2_struc(n)%clm(k)%dz(i)*250.
              clm2_struc(n)%clm(k)%h2osoi_liq(i) = 0.
           enddo
        endif
        do i = 1, nlevsoi           !soil layers
           if (clm2_struc(n)%clm(k)%t_soisno(i) <= tfrz) then       
              clm2_struc(n)%clm(k)%h2osoi_ice(i) = clm2_struc(n)%clm(k)%dz(i)* &
                   clm2_struc(n)%clm_ism*clm2_struc(n)%clm(k)%watsat(i)*denice
              clm2_struc(n)%clm(k)%h2osoi_liq(i) = 0.
              if (clm2_struc(n)%clm(k)%itypwat==istwet .or. clm2_struc(n)%clm(k)%itypwat==istice) & 
                   clm2_struc(n)%clm(k)%h2osoi_ice(i)=clm2_struc(n)%clm(k)%dz(i)*denice
           else
              if (clm2_struc(n)%clm(k)%itypwat == istsoil) then
                 clm2_struc(n)%clm(k)%h2osoi_liq(i) = clm2_struc(n)%clm(k)%dz(i)* &
                      clm2_struc(n)%clm_ism*clm2_struc(n)%clm(k)%watsat(i)*denh2o
                 clm2_struc(n)%clm(k)%h2osoi_ice(i) = 0.
                 
              elseif (clm2_struc(n)%clm(k)%itypwat==istwet .or. &
                   clm2_struc(n)%clm(k)%itypwat==istice) then 
                 clm2_struc(n)%clm(k)%h2osoi_liq(i)=clm2_struc(n)%clm(k)%dz(i)*denh2o
                 clm2_struc(n)%clm(k)%h2osoi_ice(i) = 0.                
              endif
           endif
        enddo
        
        do i = 1,nlevsoi
           if (clm2_struc(n)%clm(k)%itypwat == istsoil) then
              clm2_struc(n)%clm(k)%h2osoi_vol(i) = 0.5_r8
              clm2_struc(n)%clm(k)%h2osoi_vol(i) = clm2_struc(n)%clm(k)%h2osoi_liq(i)/&
                   (clm2_struc(n)%clm(k)%dz(i)*denh2o) &
                   + clm2_struc(n)%clm(k)%h2osoi_ice(i)/(clm2_struc(n)%clm(k)%dz(i)*denice)
           else
              clm2_struc(n)%clm(k)%h2osoi_vol(i) = 1.0_r8
           endif
           clm2_struc(n)%clm(k)%h2osoi_vol(i) = min(clm2_struc(n)%clm(k)%h2osoi_vol(i),clm2_struc(n)%clm(k)%watsat(i))
        end do
     endif
  end do

! ========================================================================
! Remaining variables are initialized by calls to ecosystem dynamics and
! albedo subroutines. 
! Note: elai, esai, frac_veg_nosno are computed in Ecosysdyn and needed
! by Fwet and SurfaceAlbedo
! Note: fwet is needed in routine clm_twostream (called by clm_surfalb)
! ========================================================================

  calday = LIS_get_curr_calday(LIS_rc)
  doalb = .true.


!  call clmlairead(n)

#if (defined DGVM)
  call iniTimeConstDGVM()
#endif

  do k = 1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     clm2_struc(n)%clm(k)%nstep = LIS_rc%nts(n)
     call EcosystemDyn (clm2_struc(n)%clm(k), doalb, .false.)
     clm2_struc(n)%clm(k)%frac_sno = clm2_struc(n)%clm(k)%snowdp/(0.1 + clm2_struc(n)%clm(k)%snowdp)  
     clm2_struc(n)%clm(k)%frac_veg_nosno = clm2_struc(n)%clm(k)%frac_veg_nosno_alb
     if(clm2_struc(n)%clm(k)%itypveg.eq.LIS_rc%bareclass.or.&
          clm2_struc(n)%clm(k)%itypveg.eq.0) then
        clm2_struc(n)%clm(k)%frac_veg_nosno = 0
     endif
     call Fwet(clm2_struc(n)%clm(k))
     call SurfaceAlbedo (clm2_struc(n)%clm(k), calday, eccen, obliqr, lambm0, mvelpp)
  end do
  return

end subroutine iniTimeVar
