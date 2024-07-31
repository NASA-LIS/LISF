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

subroutine Hydrology2 (clm)

!-----------------------------------------------------------------------
!
!  CLMCLMCLMCLMCLMCLMCLMCLMCL  A community developed and sponsored, freely   
!  L                        M  available land surface process model.  
!  M --COMMON LAND MODEL--  C  
!  C                        L  CLM WEB INFO: http://clm.gsfc.nasa.gov
!  LMCLMCLMCLMCLMCLMCLMCLMCLM  CLM ListServ/Mailing List: 
!
!-----------------------------------------------------------------------
! Purpose:
! This is the main subroutine to execute the calculation of soil/snow 
! hydrology
!
! Method:
! Calling sequence is: 
!  Hydrology2:                       surface hydrology driver
!    -> SnowWater:                   change of snow mass and snow water
!                                    onto soil
!    -> SurfaceRunoff:               surface runoff 
!    -> Infiltration:                infiltration into surface soil layer
!    -> SoilWater:                   soil water movement between layers
!          -> Tridiagonal            tridiagonal matrix solution
!    -> Drainage:                    subsurface runoff  
!    -> SnowCompaction:              compaction of snow layers
!    -> CombineSnowLayers:           combine snow layers that are thinner
!                                     than minimum
!    -> DivideSnowLayers:            subdivide snow layers that are thicker
!                                     than maximum
!    -> WetIceHydrology:             calculate hydrology for wetland and land
!                                     ice 
! Author:
!
!-----------------------------------------------------------------------
! $Id: Hydrology2.F90,v 1.6 2004/11/24 22:56:27 jim Exp $
!-----------------------------------------------------------------------

  use LIS_precisionMod
  use clm2type
  use clm2_varcon, only : denh2o, denice, istice, istwet, istsoil, spval
  use clm2_varpar, only : nlevsoi, nlevsno
  implicit none

!----Arguments----------------------------------------------------------

  type (clm1d), intent(inout)  :: clm !CLM 1-D Module

!----Local Variables----------------------------------------------------

  integer j                   ! do loop index
  real(r8) :: t_snow
  real(r8) zwice              ! the sum of ice mass of soil (kg m-2)
  real(r8) vol_liq(1:nlevsoi) ! partial volume of liquid water in layer
  real(r8) s(1:nlevsoi)       ! wetness of soil (including ice)
  real(r8) zwt                ! water table depth
  real(r8) fcov               ! fractional area with water table at surface
  real(r8) dwat(1:nlevsoi)    ! change in soil water
  real(r8) hk(1:nlevsoi)      ! hydraulic conductivity (mm h2o/s)
  real(r8) dhkdw(1:nlevsoi)   ! d(hk)/d(vol_liq)
#if (defined DGVM)
  real(r8) :: watdry
  real(r8) :: rwat     !soil water wgted by depth to maximum depth of 0.5 m
  real(r8) :: swat     !same as rwat but at saturation
  real(r8) :: rz       !thickness of soil layers contributing to rwat (m)
  real(r8) :: tsw      !volumetric soil water to 0.5 m
  real(r8) :: stsw     !volumetric soil water to 0.5 m at saturation
#endif
  real(r8) :: snowliq  
  real(r8) :: snowice
!----End Variable List--------------------------------------------------

!
! Determine the change of snow mass and the snow water onto soil
!

     call SnowWater (clm)

!
! Determine soil hydrology
!
      hk = 0.0
     if (clm%itypwat == istsoil) then
      call SurfaceRunoff  (clm, zwice, vol_liq, s, zwt, &
                           fcov)
      call Infiltration   (clm)
      call SoilWater      (clm, vol_liq, dwat, hk, dhkdw)
!      print*, 'zwice', zwice
!      print*, 'vol_liq', vol_liq
!      print*, 'zwt ',zwt
!      print*, 'fcov',fcov
!      print*, 'hk ',hk
!      print*, 'dhkdw',dhkdw
!      print*, 'dwat ',dwat
      call Drainage       (clm,  zwice, vol_liq, s,   zwt, &
                           fcov, hk,    dhkdw,   dwat )
     endif
!	print*,'check..',iam,clm%snl,nlevsno
!        clm%snl = 0
     if (clm%snl < 0) then
!        print*,'check 2..',iam, clm%snl
!
! Natural compaction and metamorphosis.
!

        call SnowCompaction (clm)

!
! Combine thin snow elements
!

        call CombineSnowLayers (clm)

!
! Divide thick snow elements
!

        call DivideSnowLayers (clm)

!
! Set empty s	print*,iam,clm%t_soisno(clm%snl+1)now layers to zero
!

!	print*,'check..',iam,clm%snl,nlevsno
        if (clm%snl > -nlevsno) then
           clm%snowage = 0.
           do j = -nlevsno+1, clm%snl
!	     print*,'h.',iam,j
              clm%h2osoi_ice(j) = 0.
              clm%h2osoi_liq(j) = 0.
              clm%t_soisno(j) = 0.
              clm%dz(j) = 0.
              clm%z(j) = 0.
              clm%zi(j-1) = 0.
           enddo
        endif

     endif
!          print*,'ps..',iam,clm%snl,clm%t_soisno
!
! Vertically average t_soisno and sum of h2osoi_liq and h2osoi_ice 
! over all snow layers for the given patch - these will be written out
! to the history tape
!
        
     if (clm%snl < 0) then
        t_snow  = 0.
        snowice = 0.
        snowliq = 0.
        do j = clm%snl+1, 0
           t_snow  = t_snow  + clm%t_soisno(j)
           snowice = snowice + clm%h2osoi_ice(j)
           snowliq = snowliq + clm%h2osoi_liq(j)
        end do
        t_snow = t_snow/abs(clm%snl)
     else
        t_snow  = spval
        snowice = spval
        snowliq = spval
     endif

!
! Update ground temperature
!

     clm%t_grnd = clm%t_soisno(clm%snl+1)
!
! Determine volumetric soil water
! 
     
     do j = 1,nlevsoi
        clm%h2osoi_vol(j) = clm%h2osoi_liq(j)/(clm%dz(j)*denh2o) &
                          + clm%h2osoi_ice(j)/(clm%dz(j)*denice)
     end do

!
! Determine ending water balance
!
     
     clm%endwb=clm%h2ocan+clm%h2osno
     do j = 1, nlevsoi
        clm%endwb = clm%endwb + clm%h2osoi_ice(j) + clm%h2osoi_liq(j)
     enddo

!     if (clm%kpatch .eq. 56063) then
!     write(*,38) clm%kpatch, clm%endwb, clm%h2ocan, clm%h2osno, clm%forc_rain, clm%forc_snow
! 38  format(i5,1x,5(f12.3,1x))
!     endif

!
! Determine wetland and land ice hydrology (must be placed here 
! since need snow updated from CombineSnowLayers)
!

     if (clm%itypwat==istwet .or. clm%itypwat==istice) call WetIceHydrology (clm)

#if (defined DGVM)
! --------------------------------------------------------------------
! Available soil water up to a depth of 0.5 m.
! Potentially available soil water (=whc) up to a depth of 0.5 m.
! Water content as fraction of whc up to a depth of 0.5 m.
! --------------------------------------------------------------------

  if (clm%itypwat == istsoil) then
     rwat = 0.
     swat = 0.
     rz = 0.
     do j = 1, nlevsoi
        if (clm%z(j)+0.5*clm%dz(j) <= 0.5) then
           watdry = clm%watsat(j) * (316230./clm%sucsat(j)) ** (-1./clm%bsw(j))
           rwat = rwat + (clm%h2osoi_vol(j)-watdry) * clm%dz(j)
           swat = swat + (clm%watsat(j)    -watdry) * clm%dz(j)
           rz = rz + clm%dz(j)
        end if
     end do
     if (rz /= 0.) then
        tsw = rwat/rz
        stsw = swat/rz
     else
        watdry = clm%watsat(1) * (316230./clm%sucsat(1)) ** (-1./clm%bsw(1))
        tsw = clm%h2osoi_vol(1) - watdry
        stsw = clm%watsat(1) - watdry
     end if
     clm%wf = tsw/stsw
  else
     clm%wf = 1.0
  end if
#endif

end subroutine Hydrology2
