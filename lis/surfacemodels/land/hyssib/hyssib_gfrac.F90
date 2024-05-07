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
!
! !ROUTINE: hyssib_gfrac.F90 
! \label{hyssib_gfrac}
!
! !REVISION HISTORY:
! 28 Apr 2002: Kristi Arsenault, Added SSiB LSM to LDAS, initial code
! 15 Feb 2004: David Mocko, Conversion from SSiB to HY-SSiB
! 28 Oct 2007: Chuck Alonge, Updates for LIS 5.0 
!
! !INTERFACE:
subroutine hyssib_gfrac()
! !USES:
  use LIS_coreMod, only : LIS_rc
  use hyssib_lsmMod 
  use hyssibveg_module

!
! !DESCRIPTION:
!  This subroutine takes vegetation greenness fraction data and the date
!  to interpolate and determine the actual value of the greenness fraction
!  for that date.  This actual value is then returned to the main program.
!  The assumption is that the data point is valid for the 16th of the given
!  month, at 00Z.
!  Gustavo will take advantage of this routine to read the other time varying
!  parameters and interpolate them to the appropriate time frame
!EOP
  implicit none
  
!=== Local Variables =====================================================
  INTEGER :: I,n,vegtyp   ! Loop counters
  INTEGER :: ITYP,IMON,ICG,IWV,ILD,IDP,IBD
!== Gustavo added for SSiB
  PARAMETER (ITYP=13,IMON=12,ICG=2,IWV=3,ILD=2,IDP=3,IBD=2)
!  REAL :: RSTPAR_v(ITYP,ICG,IWV), CHIL_v(ITYP,ICG), TOPT_v(ITYP,ICG), &
!       TLL_v(ITYP,ICG), TU_v(ITYP,ICG), DEFAC_v(ITYP,ICG),         &
!       PH1_v(ITYP,ICG), PH2_v(ITYP,ICG),  ROOTD_v(ITYP,ICG),       &
!       BEE_v(ITYP), PHSAT_v(ITYP), SATCO_v(ITYP), POROS_v(ITYP),   &
!       ZDEPTH_v(ITYP,IDP), SLOPE_v(ITYP)
!  REAL :: GREEN_v(ITYP,IMON,ICG), VCOVER_v(ITYP,IMON,ICG),           &
!       ZLT_v(ITYP,IMON,ICG), Z0_v(ITYP,IMON), DD_v(ITYP,IMON),    &
!       Z2_v(ITYP,IMON), Z1_v(ITYP,IMON), RDC_v(ITYP,IMON),        &
!       RBC_v(ITYP,IMON)

  do n=1,LIS_rc%nnest

     ! Done only once in hyssibveg_module
     !open(unit=11,file=hyssib_struc(n)%vfile,status='old', &
     !     form='unformatted')
     !read(11) rstpar_v, chil_v, topt_v, tll_v, tu_v, defac_v,   &
     !     ph1_v, ph2_v, rootd_v, bee_v, phsat_v, satco_v,   &
     !     poros_v, zdepth_v, slope_v
     !read(11) green_v, vcover_v, zlt_v, z0_v, dd_v, z2_v, z1_v,       &
     !     rdc_v, rbc_v
     !close(11)
  
!=== End Variable Definition =============================================
     do i=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
       vegtyp = hyssib_struc(n)%hyssib(i)%vegt 
       hyssib_struc(n)%hyssib(i)%vegip(1) = hyssibveg(n)%z0(vegtyp,LIS_rc%mo)
       hyssib_struc(n)%hyssib(i)%vegip(2) = hyssibveg(n)%z1(vegtyp,LIS_rc%mo)
       hyssib_struc(n)%hyssib(i)%vegip(3) = hyssibveg(n)%z2(vegtyp,LIS_rc%mo)
       hyssib_struc(n)%hyssib(i)%vegip(4) = hyssibveg(n)%dd(vegtyp,LIS_rc%mo)
       hyssib_struc(n)%hyssib(i)%vegip(5) = hyssibveg(n)%vcover(vegtyp,LIS_rc%mo,1)
       hyssib_struc(n)%hyssib(i)%vegip(6) = hyssibveg(n)%vcover(vegtyp,LIS_rc%mo,2)
       hyssib_struc(n)%hyssib(i)%vegip(7) = hyssibveg(n)%zlt(vegtyp,LIS_rc%mo,1)
       hyssib_struc(n)%hyssib(i)%vegip(8) = hyssibveg(n)%zlt(vegtyp,LIS_rc%mo,2)
       hyssib_struc(n)%hyssib(i)%vegip(9) = hyssibveg(n)%green(vegtyp,LIS_rc%mo,1)
       hyssib_struc(n)%hyssib(i)%vegip(10) = hyssibveg(n)%rbc(vegtyp,LIS_rc%mo)
       hyssib_struc(n)%hyssib(i)%vegip(11) = hyssibveg(n)%rdc(vegtyp,LIS_rc%mo)
     enddo
  enddo
end subroutine hyssib_gfrac

