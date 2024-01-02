!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!#include "LIS_misc.h"
!BOP
!
! !ROUTINE: hyssib_setvegparms
! \label{hyssib_setvegparms}
!
! !DESCRIPTION:
!  This subroutine retrieves HY-SSiB vegetation parameters
!
! !REVISION HISTORY:
! 21 Apr 2004: David Mocko, Conversion from NOAH to HY-SSiB
! 05 Sep 2007: Chuck Alonge, Updates for LIS 5.0
!
! !INTERFACE:
subroutine hyssib_setvegparms
! !USES:
  use LIS_coreMod, only : LIS_rc,LIS_surface
  use LIS_logMod, only : LIS_logunit 
  use hyssib_lsmMod
  use hyssibveg_module

! !DESCRIPTION:
!  This subroutine retrieves Hyssib veg./soil parameters. The current 
!  implementation uses a table-based lookup based on SIB classes
!  to initialize the following parameters:
!  \begin{verbatim}
!   chil    - leaf angle distribution factor
!   topt    - optimum temperature for rst calculation
!   tll     - low temp for rst calculation
!   tu      - top temp for rst calculation
!   defac   - dew factor for rst calculation
!   ph1     - stome slope factor
!   ph2     - point at which stomates close
!   rootd   - rooting depth - canopy
!   rootd   - rooting depth - ground
!   rstpar1 - PAR influence on stomatal resist. coefficients
!   rstpar2 - PAR influence on stomatal resist. coefficients
!   rstpar3 - PAR influence on stomatal resist. coefficients
!   phsat   - soil moisture potential at saturation
!   poros   - soil porosity
!   bee     - Clapp-Hornberger emperical const.
!   satco   - sat. hydraulic conductivity
!   slope   - average topographic slope %
!   zdepth1 - depth of soil layer 1
!   zdepth2 - depth of soil layer 2
!   zdepth3 - depth of soil layer 3
!  \end{verbatim}
!EOP

  implicit none

!=== Local Variables ===================================================
  integer :: i,n,vegtyp, t
!  integer :: ityp,icg,iwv,ild,idp,ibd
!  parameter(ityp=13,icg=2,iwv=3,ild=2,idp=3,ibd=2)
!  real*4 :: rstpar_v(ityp,icg,iwv),chil_v(ityp,icg),topt_v(ityp,icg),&
!       tll_v(ityp,icg),tu_v(ityp,icg),defac_v(ityp,icg),        &
!       ph1_v(ityp,icg),ph2_v(ityp,icg),rootd_v(ityp,icg),       &
!       bee_v(ityp),phsat_v(ityp),satco_v(ityp),poros_v(ityp),   &
!       zdepth_v(ityp,idp),slope_v(ityp)
!=== End Variable Definition ===========================================
  do n=1,LIS_rc%nnest 
!-----------------------------------------------------------------------
! Convert UMD Classes to SIB Classes for Each Tile
!-----------------------------------------------------------------------
     write(LIS_logunit,*) 'MSG: hyssib_setvegparms -- ',     &
          'Calling MAPVEGC to convert UMD to SIB'
     write(LIS_logunit,*) 'DBG: hyssib_setvegparms -- npatch',LIS_rc%npatch(n,LIS_rc%lsm_index)
     
     do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
        hyssib_struc(n)%hyssib(t)%vegt = &
             LIS_surface(n,LIS_rc%lsm_index)%tile(t)%vegt
        call hyssib_mapvegc(hyssib_struc(n)%hyssib(t)%vegt)
     enddo
     
!-----------------------------------------------------------------------
! Get Vegetation Parameters for HY-SSiB Model in Tile Space
! Read in the HY-SSiB Static Vegetation Parameter Files
!-----------------------------------------------------------------------
!     open(unit=11,file=hyssib_struc(n)%vfile,status='old',           &
!          form='unformatted')
!     read(11) rstpar_v, chil_v, topt_v, tll_v, tu_v, defac_v,         &
!          ph1_v, ph2_v, rootd_v, bee_v, phsat_v, satco_v,         &
!          poros_v, zdepth_v, slope_v
!     close(11)
  
!-----------------------------------------------------------------------
! Assign STATIC vegetation parameters to each tile based on the
! type of vegetation present in that tile.
!-----------------------------------------------------------------------
     do i=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
        vegtyp = hyssib_struc(n)%hyssib(i)%vegt 
        hyssib_struc(n)%hyssib(i)%vegp(1) = hyssibveg(n)%chil(vegtyp,1)
        hyssib_struc(n)%hyssib(i)%vegp(2) = hyssibveg(n)%topt(vegtyp,1)
        hyssib_struc(n)%hyssib(i)%vegp(3) = hyssibveg(n)%tll(vegtyp,1)
        hyssib_struc(n)%hyssib(i)%vegp(4) = hyssibveg(n)%tu(vegtyp,1)
        hyssib_struc(n)%hyssib(i)%vegp(5) = hyssibveg(n)%defac(vegtyp,1)
        hyssib_struc(n)%hyssib(i)%vegp(6) = hyssibveg(n)%ph1(vegtyp,1)
        hyssib_struc(n)%hyssib(i)%vegp(7) = hyssibveg(n)%ph2(vegtyp,1)
        hyssib_struc(n)%hyssib(i)%vegp(8) = hyssibveg(n)%rootd(vegtyp,1)
        hyssib_struc(n)%hyssib(i)%vegp(9) = hyssibveg(n)%rootd(vegtyp,2)
        hyssib_struc(n)%hyssib(i)%vegp(10) = hyssibveg(n)%rstpar(vegtyp,1,1)
        hyssib_struc(n)%hyssib(i)%vegp(11) = hyssibveg(n)%rstpar(vegtyp,1,2)
        hyssib_struc(n)%hyssib(i)%vegp(12) = hyssibveg(n)%rstpar(vegtyp,1,3)
        hyssib_struc(n)%hyssib(i)%vegp(13) = hyssibveg(n)%phsat(vegtyp)
        hyssib_struc(n)%hyssib(i)%vegp(14) = hyssibveg(n)%poros(vegtyp)
        hyssib_struc(n)%hyssib(i)%vegp(15) = hyssibveg(n)%bee(vegtyp)
        hyssib_struc(n)%hyssib(i)%vegp(16) = hyssibveg(n)%satco(vegtyp)
        hyssib_struc(n)%hyssib(i)%vegp(17) = hyssibveg(n)%slope(vegtyp)
        hyssib_struc(n)%hyssib(i)%vegp(18) = hyssibveg(n)%zdepth(vegtyp,1)
        hyssib_struc(n)%hyssib(i)%vegp(19) = hyssibveg(n)%zdepth(vegtyp,2)
        hyssib_struc(n)%hyssib(i)%vegp(20) = hyssibveg(n)%zdepth(vegtyp,3)
     enddo
  enddo
end subroutine hyssib_setvegparms
   
